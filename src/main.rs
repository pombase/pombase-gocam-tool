use std::{fs::File, path::PathBuf};


use clap::{Parser, Subcommand};

use pombase_gocam::parse;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    action: Action,
}

#[derive(Subcommand)]
enum Action {
    #[command(arg_required_else_help = true)]
    Stats {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    }
}

fn main() {
    let args = Args::parse();

    match args.action {
        Action::Stats { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse(&mut source).unwrap();
                println!("{}", model.id());
            }
        }
    }
}
