# Dry Lab

**Plan your PCR on your phone before you touch the wet lab.**
An offline Android app for biotechnology students and small labs: paste or open a DNA sequence and get its genes, start and stop codons, GC/AT content, PCR primers, and the full thermocycler temperature program.

Made by Paul Renaldi (Reno), B.Sc. Biotechnology, Arignar Anna College, Thiruvannamalai.

## What it does

| Screen | What you get |
|---|---|
| **Sequence** | Paste DNA or open FASTA / GenBank / text files. Length, **GC % and AT %**, base counts, DNA mass, and a GC chart that marks GC-rich spots (above 65%). |
| **Genes** | Every open reading frame on both strands: **start codon and stop codon with exact positions**, length in bp and amino acids, protein sequence and size (kDa). Genes listed inside GenBank files are shown too. |
| **Primers** | Three modes: copy a whole gene (primers on the exact start and stop), amplify around a region, or best pair anywhere. Checks Tm, GC %, 3' G/C clamp, base runs, repeats, hairpins, self-dimers, primer-dimers and extra binding sites. Optional restriction sites on primer ends for cloning. |
| **PCR** | Annealing temperature, the full program (denaturation, annealing, extension, cycles, hold) with a temperature chart, extension time from product length, touchdown option, polymerase presets, **heat advice for GC-rich and AT-rich DNA** (98 °C, DMSO, betaine, two-step PCR), gradient range, run time, and a **master mix calculator**. |
| **Tools** | Reverse complement, complement, DNA to RNA, translation in 6 frames, primer Tm, primer checker, motif search, restriction sites (18 enzymes), primer-dimer check. |

**Learn mode** explains every result in plain words, so students learn why a value matters, not just the number.

It works with **no internet** and sends nothing anywhere. Your last sequence is remembered on the phone.

## Install

1. Open the latest release in **Releases** (right side of this page) on your phone.
2. Download `DryLab-1.0.x.apk` and open it. Allow "install unknown apps" for your browser or file manager when Android asks.

Needs Android 8.0 or newer.

## How the science is done

- **Primer Tm**: nearest-neighbour thermodynamics (SantaLucia 1998 unified parameters) with salt correction and Mg²⁺ via von Ahsen 2001, the same approach as Primer3. It matches Biopython's `Tm_NN` to 0.01 °C on 900+ test primers (`tests/check_with_biopython.py`). The classroom rule 2(A+T) + 4(G+C) is shown alongside.
- **Genes**: bacterial genetic code (NCBI table 11). ATG start by default, GTG/TTG optional. Checked against Biopython translation.
- **Annealing temperature**: lower primer Tm − 4 °C for Taq, + 3 °C for Phusion/Q5-type enzymes, with the Rychlik formula as a second estimate.
- **Product Tm**: 81.5 + 16.6 log[Na⁺] + 0.41 (%GC) − 675/length.

Dry Lab is a planning and teaching tool. Always confirm results in the lab.

## Project layout

```
app/src/main/assets/www/     the app itself (HTML, CSS, JavaScript), runs offline
  js/bio.js                  biology engine: sequences, genes, Tm, primers, PCR program
  js/app.js                  screens
app/src/main/java/...        Android wrapper: file picker, share, copy
tests/test_bio.js            engine tests (node tests/test_bio.js)
tests/check_with_biopython.py   cross-check against Biopython
.github/workflows/           builds the APK and publishes a release on every push
```

Open `app/src/main/assets/www/index.html` in a browser to try it on a computer.

## Roadmap

- Tamil and other Indian-language screens
- qPCR primers and probes; agarose gel band preview
- Primer specificity search against whole genomes
- Plasmid map drawing and cloning simulation (cut, ligate, check reading frame)
- Lab notebook with photo of the gel attached to each experiment
- Shared class mode for teachers
