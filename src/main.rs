pub mod channel;
pub mod field;
pub mod merkle_tree;
pub mod polynomial;
use channel::Channel;
use chrono::Local;
use field::{Field, FieldElement};
use log::{Level, LevelFilter, Metadata, Record};
use merkle_tree::MerkleTree;
use polynomial::{interpolate_lagrange_polynomials, Polynomial};
static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;

struct ConsoleLogger;

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
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
    // 1. Generating Input
    log::info!("Generating trace");
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
    log::debug!("geneating square_trace_polynomial");
    let trace_polynomial_g = trace_polynomial.clone().compose(generator);
    let square_of_trace_polynomial_g = trace_polynomial_g.clone() * trace_polynomial_g;
    log::debug!("generating square_of_trace");
    let square_of_trace = trace_polynomial.clone() * trace_polynomial.clone();

    log::debug!("generating p_2_numerator");
    let p_2_numerator = trace_polynomial_g_2 - square_of_trace_polynomial_g - square_of_trace;
    log::info!(
        "evaluating p_2_numerator at g^1020: {:?}",
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
    // println!("{:?} : {:?} : {:?}", alpha_0, alpha_1, alpha_2);
    assert_eq!(c_p.degree(), 1023);

    log::info!("committing to composition polynomial");
    let evaluations: Vec<FieldElement> = eval_domain.into_iter().map(|h| c_p.evaluate(h)).collect();
    let merkle_tree = MerkleTree::new(&evaluations);
    // sending merkle root to channel i.e. verifier
    channel.send(merkle_tree.inner.root().unwrap().to_vec());
}

fn fibonnaci_sq(a0: FieldElement, a1: FieldElement) -> FieldElement {
    a1.pow(2) + a0.pow(2)
}

fn generate_trace(a0: FieldElement, a1: FieldElement, n: usize) -> Vec<FieldElement> {
    let mut trace = vec![a0, a1];
    for i in 2..n {
        let next = fibonnaci_sq(trace[i - 2], trace[i - 1]);
        trace.push(next);
    }
    trace
}

fn generate_eval_domain_for_trace(trace: &[FieldElement], g: FieldElement) -> Vec<FieldElement> {
    trace
        .iter()
        .enumerate()
        .map(|(i, _)| (g.pow(i as u64)))
        .collect()
}
