// in build.rs
use clap::{crate_version, load_yaml, App, AppSettings};
use clap_generate::{
    generate_to,
    generators::{Bash, Zsh},
};
use std::env;
use std::fs::File;
use std::io::prelude::*;

fn main() {
    let yaml = load_yaml!("src/cli.yaml");
    let mut app = App::from(yaml)
        .version(crate_version!())
        .setting(AppSettings::SubcommandRequiredElseHelp);

    app.set_bin_name("rustybam");
    //let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("completions/");
    let out_dir = env::var("OUT_DIR").unwrap();

    let mut file = File::create("foo.txt").expect("unable to open");
    file.write(out_dir.as_bytes()).expect("unable to write");

    generate_to::<Bash, _, _>(&mut app, "rustybam", &out_dir)
        .expect("Failed to generate bash completions");
    generate_to::<Zsh, _, _>(&mut app, "rustybam", &out_dir)
        .expect("Failed to generate zsh completions");
}
