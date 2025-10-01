#!/usr/bin/env python3

import os
from typing import List, Tuple

import matplotlib.pyplot as plt


def get_output_path(filename: str) -> str:
	"""
	Resolve an output path inside the final_statistics/EGIO_analysis directory
	relative to this script's location, creating the directory if it does not exist.
	"""
	script_dir = os.path.dirname(os.path.abspath(__file__))
	out_dir = os.path.join(script_dir, "EGIO_analysis")
	os.makedirs(out_dir, exist_ok=True)
	return os.path.join(out_dir, filename)


def prepare_data() -> Tuple[List[str], List[float], List[int]]:
	"""Return species labels, divergence times (MYA), and ortholog counts.

	Table 3.2: EGIO-identified orthologous exon pairs between human and model organisms

	Divergence times are approximate MYA values provided by the user.
	"""
	species = [
		"Mus musculus",
		"Saccharomyces cerevisiae",
		"Rattus norvegicus",
		"Drosophila melanogaster",
		"Gallus gallus",
		"Oryctolagus cuniculus",
	]

	# Estimated divergence from human (MYA)
	divergence_mya = [
		90.0,   # mouse
		1000.0, # yeast (fungi vs animals)
		90.0,   # rat
		790.0,  # fruit fly
		310.0,  # chicken (birds vs mammals)
		92.0,   # rabbit
	]

	# Orthologous exon pairs (from Table 3.2)
	ortholog_counts = [
		102_903,  # Human-Mus musculus
		0,        # Human-Saccharomyces cerevisiae
		8_311,    # Human-Rattus norvegicus
		2_329,    # Human-Drosophila melanogaster
		15_298,   # Human-Gallus gallus
		20_426,   # Human-Oryctolagus cuniculus
	]

	return species, divergence_mya, ortholog_counts


def plot_scatter(species: List[str], x_mya: List[float], y_counts: List[int], save_path: str) -> None:
	"""Create and save the scatter plot of divergence (MYA) vs ortholog counts."""
	plt.style.use("seaborn-v0_8-whitegrid")
	fig, ax = plt.subplots(figsize=(8, 5))

	ax.scatter(x_mya, y_counts, s=70, color="#1f77b4", edgecolor="black", linewidth=0.8)

	# Annotate each point with species label
	for label, x, y in zip(species, x_mya, y_counts):
		# Offset annotations slightly to avoid overlap with markers
		ax.annotate(
			label,
			(x, y),
			textcoords="offset points",
			xytext=(5, 6),
			ha="left",
			fontsize=9,
		)

	ax.set_xlabel("Evolutionary divergence from human (MYA)")
	ax.set_ylabel("Orthologous exon pairs (count)")
	ax.set_title("Evolutionary divergence vs. EGIO orthologous exon pairs")

	# Give some padding on axes
	xmin = max(0, min(x_mya) - 20)
	xmax = max(x_mya) + 100
	ymin = -max(100, max(y_counts) * 0.05)
	ymax = max(y_counts) * 1.1 if max(y_counts) > 0 else 1

	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)

	# Tight layout and save
	fig.tight_layout()
	fig.savefig(save_path, dpi=300)

	# Also show if running interactively
	try:
		plt.show()
	except Exception:
		# In headless environments, show() may fail; ignore
		pass


def main() -> None:
	species, divergence_mya, ortholog_counts = prepare_data()
	out_path = get_output_path("evolutionary_divergence_vs_orthologs.png")
	plot_scatter(species, divergence_mya, ortholog_counts, out_path)
	print(f"Saved figure to: {out_path}")


if __name__ == "__main__":
	main()
