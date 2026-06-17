import subprocess
import sys


def test_plot_cli_missing_input_returns_2(tmp_path):
    missing_csv = tmp_path / "missing.csv"
    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.plot_cli",
        "--mode",
        "continuous",
        "--pred_csv",
        str(missing_csv),
        "--title",
        "demo",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 2
    assert "predictions CSV not found" in result.stderr


def test_plot_cli_continuous_smoke(tmp_path):
    pred_csv = tmp_path / "preds.csv"
    pred_csv.write_text(
        "species,SPS,num_genes,true_phenotype\n"
        "A,0.1,10,1.0\n"
        "B,-0.2,12,2.0\n"
        "C,0.4,9,3.0\n",
        encoding="utf-8",
    )

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.plot_cli",
        "--mode",
        "continuous",
        "--pred_csv",
        str(pred_csv),
        "--title",
        "demo_plot",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "demo_plot_continuous_plot.svg").exists()
