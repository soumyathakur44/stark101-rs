pub mod field;
pub mod merkle_tree;
pub mod polynomial;
use chrono::Local;
use field::{Field, FieldElement};
use log::{Level, LevelFilter, Metadata, Record};
use merkle_tree::MerkleTree;
use polynomial::interpolate_lagrange_polynomials;
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
    let lagrange_polynomial = interpolate_lagrange_polynomials(eval_domain, trace);
    // 3. Extend
    log::info!("Extending");
    let w = generator;
    let two = FieldElement::new(2, field);
    let exp = two.pow(30) * FieldElement(3, field) / FieldElement(8192, field);
    let h = w.pow(exp.0);
    let H: Vec<FieldElement> = (0..8192).into_iter().map(|i| h.pow(i)).collect();
    let eval_domain: Vec<FieldElement> = H.into_iter().map(|h| w * h).collect();

    let evaluations: Vec<FieldElement> = eval_domain
        .into_iter()
        .map(|h| lagrange_polynomial.evaluate(h))
        .collect();
    log::info!("committing");
    // Commit to LDE
    let merkle_tree = MerkleTree::new(&evaluations);
    println!("{:?}", merkle_tree.inner.root());
    // Now we have commited to the LDE of the trace.
    // Generating constrainsts.
    // Converting polynomials to rational functions.
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
