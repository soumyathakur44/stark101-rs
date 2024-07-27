pub mod channel;
pub mod field;
pub mod merkle_tree;
pub mod polynomial;
pub mod utils;
use channel::Channel;
use chrono::Local;
use field::{Field, FieldElement};
use log::{Level, LevelFilter, Metadata, Record};
use merkle_tree::MerkleTree;
use polynomial::{interpolate_lagrange_polynomials, Polynomial};
use utils::{decommit_fri, fri_commit, generate_eval_domain_for_trace, generate_trace};
static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Debug
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!(
                "{} [{}] {}:{} - {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.module_path().unwrap(),
                record.line().unwrap(),
                record.args()
            );
        }
    }

    fn flush(&self) {}
}
// Stark 101 from stark ware written in rust.
// Fibonacci Sq Mod Prime
// a{n+2} = a{n+1}^2 + a{n}^2 mod prime
// prime = 3.2^30 + 1 = 3221225473

// Fibonacci Sq Mod 3221225473 with a{0} = 1, a{1} = x, we have a{1022} = 2338775057
fn main() {
    log::set_logger(&CONSOLE_LOGGER).unwrap();
    log::set_max_level(LevelFilter::Info);

    // Part 1: LDE and Commitment
    // LDE: Low Degree Extension in 3 steps.
    // 1. Generate Input.
    // 2. Interpolate.
    // 3. Extend i.e Evalute at many points.
    // LDE for STARK: Input a0, a1, a2, ..., a1022
    //                Evaluation domain: 1, g, g^2 g^3 ... g^1022
    //                g - element from F_p
    let prime_modulus = 3221225473;
    let field = Field::new(prime_modulus);
    let start_time = Local::now();
    // 1. Generating Input
    log::info!("Generating trace, time: {}", start_time);
    let a0 = FieldElement::new(1, field);
    let a1 = FieldElement::new(3141592, field);
    let trace = generate_trace(a0, a1, 1023);
    // get full trace paired with evaluation domain
    let generator = FieldElement::new(5, field).pow(3 * 2u64.pow(20));
    let eval_domain = generate_eval_domain_for_trace(&trace, generator);
    // 2. Interpolate
    log::info!("Interpolating");
    let trace_polynomial = interpolate_lagrange_polynomials(eval_domain, trace);
    // 3. Extend
    log::info!("Extending");
    let w = generator;
    let two = FieldElement::new(2, field);
    let exp = two.pow(30) * FieldElement(3, field) / FieldElement(8192, field);
    let h = w.pow(exp.0);

    #[allow(non_snake_case)]
    let H: Vec<FieldElement> = (0..8192).into_iter().map(|i| h.pow(i)).collect();
    let eval_domain: Vec<FieldElement> = H.into_iter().map(|h| w * h).collect();

    let evaluations: Vec<FieldElement> = eval_domain
        .clone()
        .into_iter()
        .map(|h| trace_polynomial.evaluate(h))
        .collect();

    // Commit to LDE
    log::info!("committing to LDE");
    let merkle_tree = MerkleTree::new(&evaluations);
    let mut channel = Channel::new();
    // Sending merkle root to channel i.e. verifier
    // channel acts as verifier, simulates verifier, in proof generation of non-interactive proving setup. and vice versa for proof verification
    // to convert interactive proving to non interactive proving, fiat shamir heuristic is used.
    // with fiat shamir heuristic, the channel is deterministic for both prover and verifier.
    channel.send(merkle_tree.inner.root().unwrap().to_vec());

    // Now we have commited to the LDE of the trace.
    // Generating constrainsts.
    // Converting polynomials to rational functions.
    // Build compostion polynomial and commit to composition polynomial.
    // p_0(x) = f(x) - 1 / x - g^0
    // p_1(x) = f(x) - 2338775057 / x - g^1022
    // p_2(x) = f(g^2x) - f(gx)^2 - f(x)^2 / [(x^1024 - 1)/(x-g^1021)(x-g^1022)(x-g^1023)

    // c_p = alpha_0 * p_0(x) + alpha_1 * p_1(x) + alpha_2 * p_2(x)
    log::info!("generating constraint polynomials");
    log::debug!("generating p_0");
    let p_0 = (trace_polynomial.clone()
        - Polynomial::new_from_coefficients(vec![FieldElement(1, field)]))
        / Polynomial::new_from_coefficients(vec![-generator.pow(0), FieldElement(1, field)]);
    assert_eq!(p_0.evaluate(FieldElement(2718, field)).0, 2509888982);

    log::debug!("generating p_1");
    let p_1 = (trace_polynomial.clone()
        - Polynomial::new_from_coefficients(vec![FieldElement(2338775057, field)]))
        / Polynomial::new_from_coefficients(vec![-generator.pow(1022), FieldElement(1, field)]);
    assert_eq!(p_1.evaluate(FieldElement(5772, field)).0, 232961446);

    log::debug!("generating x_1024");
    let mut x_1024 = Polynomial::new_from_coefficients(vec![FieldElement::new(0, field); 1025]);
    x_1024.coefficients[0] = -FieldElement::new(1, field);
    x_1024.coefficients[1024] = FieldElement::new(1, field);

    log::debug!("generating p_2_denominator");
    let p_2_denominator = x_1024
        / (Polynomial::new_from_coefficients(vec![-generator.pow(1021), FieldElement(1, field)])
            * Polynomial::new_from_coefficients(vec![
                -generator.pow(1022),
                FieldElement(1, field),
            ])
            * Polynomial::new_from_coefficients(vec![
                -generator.pow(1023),
                FieldElement(1, field),
            ]));

    log::debug!("generating trace_polynomial_g_2");
    let trace_polynomial_g_2 = trace_polynomial.clone().compose(generator.pow(2));
    log::debug!("generating square_trace_polynomial");
    let trace_polynomial_g = trace_polynomial.clone().compose(generator);
    let square_of_trace_polynomial_g = trace_polynomial_g.clone() * trace_polynomial_g;
    log::debug!("generating square_of_trace");
    let square_of_trace = trace_polynomial.clone() * trace_polynomial.clone();

    log::debug!("generating p_2_numerator");
    let p_2_numerator = trace_polynomial_g_2 - square_of_trace_polynomial_g - square_of_trace;
    log::debug!(
        "evaluation of p_2_numerator at g^1020: {:?}",
        p_2_numerator.evaluate(generator.pow(1020))
    );

    log::debug!("generating p_2");
    let p_2 = p_2_numerator / p_2_denominator;

    assert_eq!(p_2.degree(), 1023);
    assert_eq!(p_2.evaluate(FieldElement::new(31415, field)).0, 2090051528);

    log::info!("receiving alpha_0, alpha_1, alpha_2 from channel i.e. verifier");
    let alpha_0 = channel.receive_random_field_element(field);
    let alpha_1 = channel.receive_random_field_element(field);
    let alpha_2 = channel.receive_random_field_element(field);

    log::info!("construting compositon polynomial c_p");
    let c_p = Polynomial::new_from_coefficients(vec![alpha_0]) * p_0
        + Polynomial::new_from_coefficients(vec![alpha_1]) * p_1
        + Polynomial::new_from_coefficients(vec![alpha_2]) * p_2;

    assert_eq!(c_p.degree(), 1023);

    log::info!("committing to composition polynomial");

    // we commit to composition polynomial by evaluating it over the evaluation domain, which is the coset.
    let evaluations: Vec<FieldElement> = eval_domain
        .clone()
        .into_iter()
        .map(|h| c_p.evaluate(h))
        .collect();
    let merkle_tree = MerkleTree::new(&evaluations);
    // sending merkle root to channel i.e. verifier
    channel.send(merkle_tree.inner.root().unwrap().to_vec());

    // Now we have commited to the composition polynomial.
    // Composition polynomial will be a polynomial only if all the 3 contraints are satisfied.
    // Using F.R.I, we will show composition polynomial is close to a polynomial of low degree.
    // A function is close to a polynomial if its distance to the polynomial is small.
    // distance of a function and polynomial is measured as, f:Domain->F, Distance(f, p) => for every d in Domain, f(d) != p(d)
    // How can we prove that a composition polynomial is close to polynomial of low degree?
    // Here comes F.R.I, Fast Reed solomon Interactive oracle proofs of proximity.
    // FRI is a folding scheme, similar to what we see in FFT for breaking down the polynomial into 2 polynomials of even degree(s).
    // Prover tries to convence the verifier that, the commitment of composition polynomial is close to a low degree polynomial.
    // F.R.I Protocol:
    // 1. Receive random element `beta` from verifier.
    // 2. Apply the FRI folding scheme or FRI operator to the composition polynomial.
    // 3. Commit to the new polynomial obtained after applying FRI operator.
    // 4. Send the new commitment to the verifier.
    // 5. Repeat step 1-4, until the polynomial degree is less than accepted degree in terms of security. in this case repeat till degree is 0.
    // 6. Prover sends the result to the verifier.

    // F.R.I Operator or folding scheme:
    // from proving: function is close to a polynomial of degree < D
    // to proving: new function is close to a new polynomial of degree < D/2, where new function has half the domain size of old polynomial.
    // Example: To prove: A function is close to a polynomial of a degree < 1024, with function domain size = 8192
    //          After applying the FRI operator we need to prove the new polynomial degree < 512 with new function domain size = 4096
    // split ot even and odd powers.
    // P_0(x) = g(x^2) + x h(x^2)
    // Get random element beta from verifier
    // P_1(y) = g(y) + beta * h(y)

    // For this example, repeat steps 1-4 till degree of polynomial < 1, when domain size is 8.

    // The new evaluation domain, will be half of the old evaluation domain.
    // and new evaluation domain is first half of the old evaluation domain squared.
    // eval domain: w, w.h, w.h2, .... w.h^8091
    // new eval domain: w^2, (w.h)^2, ... (w.h^4095)^2
    // square of the first half of the old eval domian, is equal to square of second half of old eval domain. This is a cyclic group property.

    // generate fri layers and commit to the fri layers.
    log::info!("generating fri layers and fri commitments");
    let (fri_polys, _, fri_layers, fri_merkle_trees) = fri_commit(
        c_p,
        &eval_domain,
        &evaluations,
        merkle_tree.clone(),
        &mut channel,
    );

    assert!(fri_layers.len() == 11); // 11 fri layers, 8192 field elements in eval domain to 8 field elements in eval domain
    assert!(fri_layers[fri_layers.len() - 1].len() == 8); // last layer will have evaluations at 8 points
    assert!(fri_polys[fri_polys.len() - 1].degree() == 0); // and degree of last polynomial equal to zero.

    // Proof or provers work contains generating commitments and decommiting them, to convence the verifier over the integrity of the computation.
    // i)  Commitment ✅
    // ii) Decommitment -> query phase
    // Decommitment involves verifier sending random elements from evaluation domain to prover. and prover responding with decommitments to the evaluations, which involve sending merkle paths along with evaluations.
    //
    // there are a total of 12 commitments made and send over channel to verifier. The commitments made are: trace, composition polynomial, 10 fri layers.
    //
    // with each successful query and valid decommitment, verifiers confidence in the proof increases.

    log::info!("decommitting fri layers");
    decommit_fri(
        4,
        8,
        8192 - 16,
        &evaluations,
        &merkle_tree,
        &fri_layers,
        &fri_merkle_trees,
        &mut channel,
    );
    let time_taken = Local::now() - start_time;
    log::info!(
        "proof generation complete, time taken: {}ms",
        time_taken.num_milliseconds()
    );
}
