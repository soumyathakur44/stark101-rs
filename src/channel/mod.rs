use alloy::primitives::U256;

use crate::field::{Field, FieldElement};
pub struct Channel {
    pub state: String,
    pub proof: Vec<Vec<u8>>,
}

impl Channel {
    pub fn new() -> Channel {
        Channel {
            state: String::new(),
            proof: Vec::new(),
        }
    }

    pub fn send(&mut self, s: Vec<u8>) {
        log::debug!("sending to channel");
        let data_for_digest = format!("{}:{:?}", self.state, s.clone());
        self.state = sha256::digest(data_for_digest);
        // in stark101 from starkware, we push parent function and s into the proof.
        // there is no straight forward way to know parent function in rust.
        self.proof.push(s)
    }

    pub fn receive_random_field_element(&mut self, field: Field) -> FieldElement {
        let received_int = self.receive_random_int(0, field.0 - 1, true);
        log::debug!("received_int: {}", received_int);
        FieldElement::new(received_int, field)
    }

    pub fn receive_random_int(&mut self, min: u64, max: u64, show_in_proof: bool) -> u64 {
        // convert state to hexadecimal number
        let num = (U256::from(min) + U256::from_str_radix(&self.state, 16).unwrap())
            % U256::from(max - min + 1);
        let state = self.state.clone() + num.to_string().as_str();
        self.state = sha256::digest(state);
        if show_in_proof {
            self.proof.push(num.to_string().into_bytes());
        }
        num.into_limbs()[0]
    }
}
