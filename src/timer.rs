//! Simple timing utility

use std::time::Instant;

pub struct Timer {
    name: &'static str,
    start: Instant,
}

impl Timer {
    pub fn new(name: &'static str) -> Self {
        #[cfg(feature = "profile")]
        eprintln!("  * {}", name);
        
        Timer {
            name,
            start: Instant::now(),
        }
    }

    pub fn stop(&self) {
        let duration = self.start.elapsed();
        #[cfg(feature = "profile")]
        eprintln!("  * {} {:?}", self.name, duration);
        let _ = duration; // Suppress unused warning when not profiling
    }

    pub fn print(msg: &str) {
        #[cfg(feature = "profile")]
        eprintln!("  * {}", msg);
        let _ = msg; // Suppress unused warning when not profiling
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        // Automatically print timing on drop if not manually stopped
    }
}
