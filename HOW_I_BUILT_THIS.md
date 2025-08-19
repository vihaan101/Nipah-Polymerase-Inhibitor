# Scientific Implementation Journey: A Retrospective Analysis

**Author:** Vihaan Agrawal
**Subject:** Methodological Evolution and Engineering Decisions
**Date:** January 2026

---

## 1. Introduction

This isn't the clean version of the project where everything worked on the first try. This is the log of the engineering hurdles I hit, the specific failures that almost killed the project, and the reasoning I used to fix them. My goal was simple: design a Nipah Virus inhibitor that works even against the W730A resistance mutation. The path to get there was... messy.

---

## 2. Phase 1: Establishing the Ground Truth

### The Problem of Structure
I started out trying to use **Homology Modeling**—basically guessing the Nipah structure by looking at its cousin, the Measles virus. It seemed like a good shortcut, but I realized pretty quickly that "guessing" wasn't good enough. Homology models often have uncertainties of >2.0 Å in side-chain positioning. In an allosteric pocket, being off by 2 Å is like trying to unlock a door with the wrong key.

So I scrapped the models and switched to the newly released **Cryo-EM data (PDB: 9KNZ)**. Since this is experimental data derived from actual electron density maps, it gave me a "ground truth" coordinate system to start from.

### The "OpenBabel Truncation" Bug
Once I had the structure, I needed to convert it to PDBQT format for the physics engine. I threw it into the standard MGLTools pipeline and moved on. Days later, I noticed my docking results were weird. I opened the file and realized the entire C-terminal domain (residues 1850+) was gone.

It turned out the massive L-Protein (~2000 residues) was triggering a silent buffer overflow in the standard libraries. The tool didn't crash; it just stopped writing.

I had to write a custom script (`splice_receptor.py`) to split the protein into "N-Term" and "C-Term" chunks, convert them separately to bypass the buffer limit, and then surgically stitch them back together.

### The Crystal Water Trap
I also messed up the water handling at first. The crystal structure came with a bunch of $H_2O$ molecules trapped inside, and I assumed they were important structural bridges. Big mistake. My early docking runs were terrible—the drugs refused to bind because I was essentially asking them to squeeze into a pocket filled with concrete.

Then I realized I was ignoring **entropy**. In a real cell, those waters aren't frozen; they're chaotic. When a drug enters, it kicks them out (desolvation), which actually *helps* binding. Once I deleted the explicit waters, the drugs finally fit.

---

## 3. Phase 2: Designing the Physics Engine

### The "Vacuum Hole" Fallacy
I wanted to simulate **Induced Fit** (the protein hugging the drug), so I turned on "Flexible Receptor Docking". It sounded smarter. It wasn't.

When I mutated the large Tryptophan (W730) to a tiny Alanine, it left a physical hole in the protein. The flexible docking algorithm saw free real estate and twisted the adjacent side-chains *into* this void to touch the drug. I was getting amazing affinity scores, but they were lies.

Real proteins don't collapse like that on a nanosecond timescale; the beta-sheet backbone holds them rigid. I was simulating a physics impossibility. I had to pivot to a **Rigid Backbone Assumption**, forcing the drug to find a fit without "cheating" by relying on the protein to deform into empty space.

### Scoring Functions & Search Depth
I also ditched the default Vina scoring function. It over-prioritizes hydrogen bonds, but my target pocket is deeply hydrophobic (oily). I switched to **Vinardo**, which has a better dispersion term for lipophilic interactions.

Docking is basically rolling dice (Monte Carlo search). With the default settings (exhaustiveness=8), I was seeing my scores jump around by ~0.5 kcal/mol. That's too much noise when you're looking for subtle resistance effects. I cranked the exhaustiveness up to **16**, effectively doubling the search depth. It slowed everything down, but the variance dropped to near zero.

---

## 4. Phase 3: Candidate Selection

### The High Affinity Trap
Early on, I found a molecule called **ERDRP-0519**. It had a massive binding affinity (-6.69 kcal/mol). I almost declared victory right there.

Then I checked the literature: ERDRP-0519 is known to fail against the W730A mutant.

That was a wake-up call. **Affinity $\neq$ Resistance**. A drug can stick like glue to the Wild-Type but fall off the Mutant. I stopped looking for the "strongest" number and started calculating the "Delta" ($\Delta \Delta G$). I prioritized drugs where the score *didn't change* between the two states, even if the absolute score was slightly lower. Reliability beats peak performance.

### The Geometric "Ghost Clash"
Even with the Delta metric, I was getting false positives. Some drugs were scoring well against the Mutant by docking exactly where the deleted Tryptophan used to be.

These "Mutant-Selective" drugs are useless because they'd crash into the Tryptophan in the Wild-Type virus. To fix this, I wrote a geometric filter using `NumPy`. I calculated the distance between the drug and the "Ghost Atoms" of the missing Tryptophan. If a drug got within 1.5 Å of the ghost, I killed it. This forced the drug to respect the volume of *both* states.

---

## 5. Phase 4: Validating Real-World Use

Finally, I had to look at my winner, **BMS-986205**. It violated Lipinski's Rule of 5 (LogP was 6.58). In a pure chem-informatics class, this would be a "fail".

But context matters. BMS-986205 isn't a theoretical molecule; it's physically in Phase 3 clinical trials for cancer. That means real human beings are swallowing it and it's working. The empirical fact that it has a formulated delivery mechanism overrides the theoretical rule of thumb. I decided to keep it, citing the clinical data as proof of bioavailability.
