import run_differential_analyses
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('het_score_path', type=str, help='Please specify path to heterogeneity score.')
    return parser
        

if __name__ == '__main__':
    args = _get_parser().parse_args()
    run_differential_analyses.run(args.het_score_path)


















