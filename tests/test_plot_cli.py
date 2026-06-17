import subprocess
import sys


def test_plot_cli_binary_forwards_custom_phenotype_names(tmp_path, monkeypatch):
    from esl_psc_cli import plot_cli

    pred_csv = tmp_path / "preds.csv"
    pred_csv.write_text(
        "species,SPS,num_genes,input_RMSE,true_phenotype\n"
        "A,0.5,10,0.1,1\n"
        "B,-0.5,10,0.2,-1\n",
        encoding="utf-8",
    )
    seen = {}

    def fake_rmse_range_pred_plots(pred_csv_path, title, pheno_names=None, min_genes=0, plot_type="violin"):
        seen["pred_csv_path"] = pred_csv_path
        seen["title"] = title
        seen["pheno_names"] = pheno_names
        seen["min_genes"] = min_genes
        seen["plot_type"] = plot_type

    monkeypatch.setattr(plot_cli.ecf, "rmse_range_pred_plots", fake_rmse_range_pred_plots)

    rc = plot_cli.main([
        "--mode", "violin",
        "--pred_csv", str(pred_csv),
        "--title", "demo_plot",
        "--pheno_name1", "C4",
        "--pheno_name2", "C3",
    ])

    assert rc == 0
    assert seen["pheno_names"] == ("C4", "C3")
    assert seen["plot_type"] == "violin"


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
