"""Stress-test the storage pipeline under harsher error rates than the
defaults, to see where each codec's redundancy budget breaks down.

Run: python examples/stress_test.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore import DNAStorage


def trial(codec_name: str, parity_fraction: float, fountain_redundancy: float,
          substitution_rate: float, dropout_rate: float, trials: int = 20) -> float:
    successes = 0
    data = os.urandom(2000)
    for _ in range(trials):
        dna = DNAStorage(
            codec=codec_name,
            parity_fraction=parity_fraction,
            fountain_redundancy=fountain_redundancy,
            synthesis_kwargs={"substitution_rate": substitution_rate, "dropout_rate": dropout_rate},
            sequencing_kwargs={"substitution_rate": substitution_rate, "read_dropout_rate": dropout_rate},
        )
        dna.store(data, name="f.bin")
        try:
            successes += dna.retrieve("f.bin") == data
        except Exception:
            pass
    return successes / trials


def main():
    print("Stress-testing at elevated error rates (defaults are ~0.02% subst, 1% dropout)\n")
    scenarios = [
        ("mild (defaults)", 0.0002, 0.01),
        ("moderate", 0.001, 0.03),
        ("harsh", 0.003, 0.08),
    ]

    print(f"{'scenario':<20} {'naive (0.4 parity)':<20} {'rotating (0.4 parity)':<22} {'fountain (2.5x)':<16}")
    for label, subst, dropout in scenarios:
        naive_rate = trial("naive", 0.4, 2.5, subst, dropout)
        rotating_rate = trial("rotating", 0.4, 2.5, subst, dropout)
        fountain_rate = trial("fountain", 0.4, 2.5, subst, dropout)
        print(f"{label:<20} {naive_rate*100:<20.0f} {rotating_rate*100:<22.0f} {fountain_rate*100:<16.0f}")

    print("\nValues are % of trials that round-tripped successfully.")
    print("If a scenario drops to 0%, that redundancy budget is insufficient --")
    print("increase parity_fraction (naive/rotating) or fountain_redundancy (fountain).")


if __name__ == "__main__":
    main()
