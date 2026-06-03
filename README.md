# Central Limit Theorem for Means

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

The Central Limit Theorem (CLT) says that if you take the mean of a large
enough random sample, that mean will be approximately normally distributed —
even if the original data are not. This app lets you see that happen.

Choose a parent distribution (normal, uniform, right-skewed, or left-skewed),
set the sample size and number of samples, and step through three tabs:

1. **Population Distribution** — the shape of the population you are sampling from
2. **Samples** — eight individual samples drawn from that population
3. **Sampling Distribution** — the distribution of means across all samples,
   with an ideal normal curve overlaid for comparison

## How to use

1. Choose a parent distribution and adjust its parameters.
2. Set the sample size (n) — try small values like 2 or 5, then increase to 30+.
3. Set the number of samples (k) — more samples give a smoother sampling distribution.
4. Step through the three tabs to follow the logic of the CLT.

## Learning goals

- Understand what a sampling distribution is and how it is built
- See that means of repeated samples are always more normally distributed
  than the raw data
- Observe how sample size controls how quickly the CLT kicks in

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).