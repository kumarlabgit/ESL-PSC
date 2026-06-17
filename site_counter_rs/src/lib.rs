#[allow(dead_code)]
#[path = "main.rs"]
mod backend_impl;

pub use backend_impl::run_backend_stdio;
