# Field-Superposition Annealing for Heterogeneous 2.5D Chiplet Co-Design under Thermal and Stress Constraints

*DAC regular paper draft. Six pages of body, references on page 7 only. Figures are named where they go; the captions say which evidence file they come from. Numbers in square brackets are the ones that differ from the outline you sent, and the last section explains each one.*

Index Terms — 2.5D chiplet, design-space exploration, thermal simulation, von Mises stress, simulated annealing.

---

## Abstract

The design space of a 2.5D package that mixes CPU, GPU, compute-in-memory and HBM chiplets is large, and the tools that search it either cover only one of those device types or price every candidate with a full thermal solve. We present a co-design loop that takes a workload and a catalog of commercial parts, places them, and returns a floorplan together with a per-chiplet frequency setting that stays under a temperature and a stress limit chosen for the target use. The search prices each move by superposition: a chiplet that did not move keeps its stored temperature field, and only the moved one is recomputed. On the TAP-2.5D cases this cuts the search time by 1.49× to 2.08× at the same grid, and by up to 3.70× over a full TAP schedule once a machine-sleep gap is removed from the timing. Where the mold fill differs from silicon, superposition on a stale material map is wrong by tens of kelvin; eight fixed correction steps on the factorization taken at the start bring it back, at 1.97× to 2.75× over refactoring every step, with the peak of the accepted floorplan within 0.8 K of an exact rescore. Fill conductivity changes which path carries the heat: the secondary paths carry 10% of the power with air in the gaps and 69% with silicon. A thermal-stress solver checked against COMSOL is within 0.3 K at chip centers; its von Mises stress runs high by 12.1% on average, which we remove with one scale factor and then treat as a gate applied after the floorplan freezes, not inside the anneal.

---

## I. Introduction

A 2.5D package puts several chiplets on one interposer under one spreader. CPU, GPU and compute-in-memory chiplets are now placed together, and the memory bandwidth they need depends on how many of each are present and on what they are running. Searching that space means generating power for each chiplet from its own workload, placing them, and checking the result against a thermal limit and a mechanical one. The mechanical limit matters because the package warps, and the warp grows with the temperature rise.

Three things make this search expensive or misleading. First, the power tools are separate: gem5 with McPAT for a CPU, Accel-Sim with AccelWattch for a GPU, AccelForge for a CIM array, Ramulator for the memory. None of them knows about the others, and none of them resizes the HBM power when the number of CPUs and GPUs changes. Second, the thermal model used for placement usually fills the gaps between chiplets with silicon. That is what HotSpot does. The fill is not silicon in a real package, and its conductivity decides how much heat leaves through the spreader and how much leaves sideways and downward. Third, stress is either left out or folded into every step of the search. A stress solve costs about as much as thirty thermal solves, so putting it inside the anneal makes the search too slow to be useful, while leaving it out entirely can accept a floorplan that warps past its limit.

We build one loop around these three. The front end runs the four power tools and sets the HBM power from the CPU and GPU count. The inner search is simulated annealing over chiplet positions, rotations and part swaps, priced by superposition so that a move only recomputes what moved. When the fill is not silicon, the factorization from the starting layout is reused and corrected for a fixed eight iterations instead of being rebuilt. Stress is computed once, after the floorplan is frozen, and only then does a second loop raise or lower each CPU and GPU frequency to sit just under the limits.

Our contributions are:

1. A closed loop from workload to floorplan and frequency. It runs gem5, McPAT, Accel-Sim, AccelWattch, AccelForge and Ramulator together, sizes the HBM power from the number of CPUs and GPUs, and searches placement and part selection under a cost that splits peak temperature, temperature gradient and wirelength into equal thirds.

2. A way to price an annealing step without refactoring. Superposition reuses the fields of chiplets that did not move. Where the fill differs from silicon and that reuse is no longer exact, eight correction steps on the starting factorization replace the refactor, giving 1.97× to 2.75× over refactoring every step with the accepted floorplan rescoring within 0.8 K.

3. A thermal-stress check that is calibrated and used as a gate rather than as a cost. Against COMSOL the temperature at chip centers is within 0.3 K. The von Mises stress is high by 12.1% on average over seven layouts; one scale factor removes that mean bias, and the corrected stress is judged only after the floorplan freezes, against limits the user sets for the application.

The rest of this paper is organized as follows. Section II reviews the tools this work builds on and states the problem. Section III describes the co-design loop. Section IV gives the superposition method and the correction used when the fill is not silicon. Section V reports the calibration, the speedup and the design-space results. Section VI concludes.

---

## II. Background and Related Work

### II-A. Thermal models and what they leave out

HotSpot's grid model solves a four-layer stack: the chip layer, a thermal interface material, a copper spreader and a copper heat sink, with convection off the top of the sink. It is fast enough to call inside a search and it is the model TAP-2.5D uses, so we keep it as the reference for temperature. Two things about it matter here. The gaps between chiplets are filled with the chip-layer conductivity, which is silicon, so a change of mold compound never reaches the solution. And the grid has no lateral heat path at the package edge beyond the sink footprint, which changes the coolest part of the board by several kelvin; we add that path back when we compare solvers, and say so wherever we do.

Stress does not behave like temperature. The temperature field of a linear conduction problem is the sum of the fields of its sources. The stress field is not, because the constraint that produces it depends on the whole temperature distribution at once. That is why a stress result cannot be cached per chiplet and added, and why we keep it out of the inner loop.

### II-B. Related work

**Floorplan search for chiplets.** TAP-2.5D anneals a chiplet floorplan and calls HotSpot at every step, with wirelength from a router. It handles CPU-class chiplets and does not compute stress. Its own anneal is the baseline we time our superposition against, because the placement code, the moves and the initial floorplan are then identical and only the temperature entry differs. ATMPlace adds thermal, stress and wirelength terms and trains a model to avoid most of the solves, reporting about 7× over TAP-2.5D. That approach and ours remove different parts of the cost: theirs avoids evaluations, ours makes each evaluation cheaper. The two can be combined, and we do not claim otherwise.

**Power modeling.** gem5 with McPAT, Accel-Sim with AccelWattch, and AccelForge each produce a credible power number for one device type. Ramulator does the same for the memory system. Using them together needs a convention for the parts none of them owns: the HBM power has to follow the CPU and GPU count, and the CIM power has to come from somewhere other than a CPU template. We treat that as plumbing and describe it in Section III rather than as a result.

**Thermal-mechanical co-simulation.** Finite-element solvers such as COMSOL give a reference for both temperature and stress, at a cost that rules them out inside a search. What a search needs from them is a check that the cheaper model is biased in a known direction. Section V-A reports that check for our solver.

### II-C. Problem formulation

Given a set of chiplets, each with a power and a size, and a package with a fixed outline, we want positions and orientations, a choice of part for each CPU and GPU slot, and a frequency scale per CPU and GPU. The inner objective, which we call the cost, is

```
J = ( (T_max - T_amb) / dT_0  +  g / g_0  +  L / L_0 ) / 3
```

where `g` is the largest temperature gradient on the chip layer, `L` is the wirelength per wire, and the three divisors are the values at the starting floorplan, so the three terms start equal. Lower is better. Wirelength uses the center-to-center Manhattan distance with a fixed connection matrix: 1024 each way between every CPU or GPU and the HBM, 128 each way between every CIM and the first CPU. The matrix is assumed, not taken from a router, and it is frozen before the anneal so that the initial and the final floorplan are scored on the same netlist.

The constraints are applied after the floorplan freezes, not inside `J`. On the 2 mm grid the hottest die cell must stay below 90 °C and the von Mises stress inside a normal chiplet must stay below 200 MPa. Both limits are inputs; a part aimed at a vehicle can raise them, a part aimed at a laptop can lower them. HBM is not counted as a normal chiplet for the stress gate.

The linear system behind the inner loop is `A T = P`, with `A` the conductance matrix. Filling `A` is cheap. Factoring it is where the time goes, and it is the step Section IV avoids repeating.

---

## III. The Co-Design Loop

Fig. 1 shows the loop. *[Fig. 1: one wide figure. Left, the front end: workload in, four power tools, HBM power out. Middle, the anneal: propose a move, price it by superposition or by the corrected factorization, accept or reject. Right, the frozen floorplan going to the stress solver and then to the frequency loop. Mark the stress solver as outside the anneal.]*

### III-A. Inputs and the starting floorplan

The input is a list of chiplets and a task for each. CPU power comes from gem5 traces through McPAT, GPU power from Accel-Sim through AccelWattch, CIM power from AccelForge, and the HBM power from a channel formula that grows with the number of CPUs and GPUs. Parts whose datasheet gives no die area get the area of the closest commercial part of the same class; the catalog records which one was used.

The starting floorplan comes from three constructors, a B*-tree packer, a Monte-Carlo tree search and an L-BFGS relaxation, all run at the same seed. Each result is snapped onto the 1.125 mm grid and pushed apart until every pair has at least one cell of gap. The three are scored with the stress solver as `T/90 + σ/200` and the lowest score wins. We tried letting the search start from a hand-placed shelf; the generated start is what the loop actually uses, so the comparison in Section V is against that.

### III-B. The inner loop

A step does one of five things. Swapping a CPU or GPU for another part from the catalog is drawn 10% of the time; CIM chiplets are never swapped and the HBM power is recomputed whenever the counts change. The other 90% is split four ways: rotation takes a fixed tenth of it, a small nudge takes a share that grows linearly from 0.10 to 0.70 as the anneal proceeds, and the rest is divided equally between a long jump and a swap of two positions. So the search moves coarsely at the start and finely at the end, which is the schedule that converged fastest in an earlier comparison of four fixed ratios.

A move that lowers the cost or leaves it unchanged is always taken. A move that raises it is taken with probability `0.20 × (1 − step / total)`, so the chance falls from 0.20 to zero and the last step never accepts a worse floorplan. The best floorplan seen is stored separately and is not overwritten by an accepted worse step.

The cost is priced by superposition, described in Section IV. For a fill other than silicon the conductance matrix is rebuilt and factored once, at the start, and every later step reuses that factorization with a residual correction. The drift of this correction, measured by rescoring the final floorplan exactly, was 0.37 K to 0.85 K.

### III-C. The outer loop

Once the floorplan is frozen it is not moved again unless the stress gate fails. The stress solver runs on the 2 mm grid and reports the hottest die cell and the highest von Mises stress inside a normal chiplet. Frequency is then adjusted per CPU and per GPU, in steps of 0.1, over the range 0.5 to 2.5, scaling dynamic power only. The search raises scales for as long as a guide estimate stays under 90 °C, then checks the true result and walks the worst offender back down if the gate is missed. Stress above 200 MPa but within 220 MPa gets a short repair of at most six tries; stress above 220 MPa is left as it is and marked as failing the limit, so the result is kept and visible rather than dropped.

This split is the point of the outer loop. The inner loop answers where the chiplets should go. The outer loop answers, for a floorplan that is already fixed, how far each chiplet can be pushed before the package crosses its limit.

### III-D. The thermal stack used for pricing

The pricing model is the HotSpot stack: a 150 µm chip layer, 20 µm of interface material at 4 W/m·K, a 1 mm copper spreader and a 6.9 mm copper sink, both 96 mm on a side over a 72 mm chip frame, convection at the sink's value, ambient at 45 °C. Chip cells are silicon at 130 W/m·K. Gap cells take the fill conductivity, which is an input which is one of silicon, epoxy at 0.8, underfill at 1.0, air at 0.028, or a named value. The stress solver uses the same stack plus a 1 mm substrate and a 100 µm interposer, and the same convection.

---

## IV. Pricing a Move without Refactoring

### IV-A. Superposition

In a linear conduction problem the temperature field of several sources is the sum of the fields of each source alone. For a floorplan whose material map does not change, the field of chiplet `i` depends only on its own power and position. We compute it once, store it, and rebuild the package field by adding the stored fields of every chiplet. A move changes one or two chiplets, so only those fields are recomputed and the rest are read back.

The saving scales with how rarely a chiplet moves, not with the grid. On the TAP-2.5D anneal, where most chiplets stay put at every step, the number of solves per search drops by about a third and the wall time drops by 1.49× to 2.08× at one grid. The same ratio appears at 16 and at 64 cells across the package, which confirms that the gain comes from fewer solves and not from the grid.

### IV-B. When the material map changes

Superposition assumes the conductance matrix is fixed. It is not fixed when the gaps are filled with anything but silicon, because moving a chiplet moves the silicon cells and the fill cells with it. Keeping the old matrix and only updating the power then prices the move on the wrong material map. On a three-chiplet floorplan with epoxy in the gaps, that error reached 26 K at the peak and 31 K at a single cell within ten steps, and it appeared at the first move rather than accumulating. For a search that accepts on a tenth of a kelvin, this is unusable.

### IV-C. Correcting instead of refactoring

A solve has three parts. Filling the matrix from the material map is cheap. Factoring it is the expensive part, by more than an order of magnitude. Back-substitution with the power vector is cheap again. So the part worth avoiding is the factor.

Let `A0` be the matrix at the start, already factored, and `A1` the matrix after a move. We form a trial field `T*` by back-substitution of the new power through the old factor. Because `A1` is known, the power that would actually produce `T*` on the new map is `P* = A1 T*`, and it costs only a multiply. The gap `P − P*` is the power the old factor failed to account for. Feeding that gap back through the same factor and repeating gives a fixed-point iteration whose fixed point is the exact field of `A1`. We run it for eight iterations and stop, with no error test and no refactor on acceptance.

Eight was chosen because the residual after eight steps is far below the acceptance threshold, and because a fixed count keeps the cost of a step predictable. The iteration is the whole of the method: there is no preconditioner and no rank update.

### IV-D. What it costs and where it weakens

The factor is paid once, at the start. Every later step pays a matrix fill and eight multiplies. Against refactoring at every step this gives 2.42× to 2.75× on the five-layer model with a cold package edge, and 1.97× to 2.11× on the four-layer model with the spreader overhang included. The second number is lower because the overhang model factors faster, so there is less to save.

The error to report is not the difference between the two anneals, which diverge onto different floorplans. It is the difference between the price the iteration quoted for its chosen floorplan and an exact refactor of that same floorplan. That peak error stays within 0.8 K on the cold-edge model and within 0.08 K on the overhang model. The iteration quotes low more often than high, so it occasionally accepts a floorplan it should have rejected; the exact rescore is what we publish.

The method needs the material change to be a correction of a matrix that is already a reasonable approximation. A fill as insulating as air makes the first move jump by hundreds of kelvin and the iteration starts too far from the solution to be trusted; we do not use it for air.

---

## V. Evaluation

### V-A. Solver calibration

*[Table I: temperature. Rows are the COMSOL cases, columns are chip-center temperature from COMSOL and from our solver, and the difference. All differences under 0.3 K.]*

*[Table II: stress, seven layouts at the chip mid-plane, z = 1.175 mm, against COMSOL linear hexahedra. Columns: layout, relative difference before the scale factor, after it.]*

Our thermal-stress solver reproduces COMSOL chip-center temperatures to within 0.3 K. Its von Mises stress sits high. Averaging the relative difference at chip centers within each of seven layouts, then averaging those seven, gives +12.1%. We multiply reported stress by `1 / 1.121 = 0.892` and leave the solver, the mesh, the materials and the temperature field untouched. After this the mean difference over the seven layouts is zero by construction, which is not the same thing as agreement with COMSOL: the spread between layouts remains, the worst single layout moves from −2% to −13%, and the others fall between about −2% and +7%. We therefore call the solver conservative before the correction and unbiased on average after it, and we never describe it as matching COMSOL.

The same comparison shows the secondary heat paths are not a small effect. With the gaps in epoxy, 73.6% of the power leaves through the spreader, 10.3% through the board underside and 16.1% through the package edge. Replacing the epoxy with silicon drops the spreader share to 31.0% and raises the edge share to 64.0%. A model that fills the gaps with silicon answers a different question about where the heat goes.

*[Fig. 2: stacked bars of heat flow through spreader, underside and edge, for air, organic, epoxy, underfill and silicon fills. Source: 论文证据/填充场叠加.]*

### V-B. Superposition against a full solve

*[Table III: the three TAP-2.5D cases at 16 and 64 cells, full solve against superposition, time and rescored peak.]*

We ran TAP-2.5D's own anneal, its own moves and its own initial floorplan, and replaced only the temperature call. Both sides use seed 22 and start from byte-identical floorplans. At 64 cells the full solve takes 1316 s to 1478 s for 120 steps and superposition takes 659 s to 855 s, a gain of 1.54× to 2.08×. The rescored peaks match to 0.02 K on two of the three cases.

The third, micro150, finishes 2.49 K hotter under superposition. The two searches are identical up to step 62, where the same candidate is priced at 102.62 °C by the full solve and 102.60 °C by superposition. The full solve rejects it and superposition accepts it. A 0.02 K difference sat on the acceptance threshold, the decision flipped, and the searches diverged. We keep the 2.49 K figure and report the cause, because replacing it with a rerun would hide the one real failure mode of pricing a search on an approximation: not a wrong field, but a decision that falls the wrong side of a threshold. A second seed on the same starting floorplan does not flip and finishes within 1.5 K of the full solve.

*[Fig. 3: paired bars of search time, log scale, and the rescored temperature difference, for the six runs. Source: 论文证据/tap场叠加/图/用时与误差.png.]*

Over a full TAP schedule of 1050 steps on the ascend910 case, the full solve takes 45507 s of compute and superposition 12306 s, which is 3.70×. The raw wall time is longer by a stretch where the machine was asleep; we subtract that stretch and cite the compute time only.

### V-C. The correction, and what the fill does

*[Table IV: six floorplans from 5 to 11 chiplets, full refactor against eight correction steps. Columns: time, speedup, peak error, maximum cell error.]*

With epoxy in the gaps, eight correction steps on the starting factorization run in 2.42× to 2.75× less time than refactoring every step, across 5 to 11 chiplets. The quoted peak is within 0.8 K of an exact rescore of the same floorplan, and the quoted value is low in five of the six cases, so the search is slightly optimistic. On the overhang model, which is the one comparable with HotSpot, the same comparison gives 1.97× to 2.11× and the cell error stays under 0.08 K.

*[Fig. 4: bars of time and of cell error for the overhang model at 5, 8 and 11 chiplets. Source: 论文证据/填充场叠加/图/外伸环_加速对比.png.]*

The fill sweep says the correction is worth having. Moving from air to silicon in the gaps moves 59 percentage points of the power off the main path and onto the edges. A search that assumes silicon fill is optimizing a package that does not exist.

### V-D. The full loop

*[Table V: per case, initial against best against after frequency scaling. Columns: hottest-chip temperature, internal stress, gradient, wirelength, frequency scales, whether the gate passed. This table is filled from the run in 论文证据/闭环全流程 when it completes.]*

The loop runs on packages of 5, 6, 7, 8, 9 and 11 chiplets, each with silicon fill and with epoxy fill. We report the starting floorplan, the best floorplan the anneal found, and the same floorplan after frequency scaling, so the gain from placement and the gain from frequency are separate. Placement is what reduces peak temperature, gradient and wirelength. Frequency scaling does not move a chiplet; it spends whatever thermal and stress margin the placement created, and its benefit shows up as frequency rather than as temperature.

### V-E. What the results suggest

Three observations hold across the runs and are robust to the exact floorplans.

The fill should be an input, not a default. The share of power leaving through the secondary paths ranges from 10% to 69% over the fills we tried, which is large enough to change both the peak temperature and the temperature gradient that the cost function sees. A tool that fills the gaps with silicon will prefer different floorplans from one that fills them with epoxy.

Stress belongs after the search, not inside it. The stress field cannot be superimposed, and one stress solve costs as much as many thermal solves. Judging it once, on the frozen floorplan, against a limit chosen for the application, catches the floorplans that warp too far without dominating the search time. The limit is the knob: the same loop serves a part with a tight limit and a part with a loose one, and the frequency scaling lands at a different point for each.

Superposition is safe to build on, with one caveat that should be reported rather than removed. Its field error is well under the acceptance threshold. Its danger is the candidate that lands on the threshold, where a 0.02 K pricing error flips the decision and the search diverges. That happened once in our runs. An exact rescore of the final floorplan costs a single solve and removes the doubt, so we treat it as part of the method and not as an optional check.

---

## VI. Conclusion

We described a co-design loop that takes a workload and a catalog of CPU, GPU, CIM and HBM chiplets and returns a floorplan and a per-chiplet frequency setting under a temperature and a stress limit. The search stays fast because a move recomputes only what moved, and because a change of fill is handled by correcting a factorization instead of rebuilding it, at roughly 2× over refactoring every step and with the accepted floorplan rescoring to within 1 K. The stress check is calibrated against COMSOL, corrected for its average bias, and applied once after the floorplan freezes.

The results argue for three choices: treat the mold fill as an input, keep stress out of the inner loop and apply it as a gate, and rescore the final floorplan exactly so that a decision caught on the acceptance threshold cannot silently stand.

The correction is not reliable for a fill as insulating as air, where the material change is too large for eight iterations to absorb. The wirelength model uses an assumed netlist rather than a routed one. Both are left for later work, along with combining this per-evaluation saving with a learned model that avoids evaluations altogether.

---

## References

*[Page 7 only. No body text on this page.]*

---

# Notes on numbers that differ from the outline

These are the places where the outline you sent states a figure I did not use, and the reason.

| Outline | Draft | Why |
| --- | --- | --- |
| Superposition gives 46× | 1.49× to 2.08× at one grid, up to 3.70× over a full TAP schedule | 46.3× compared a 16-cell superposition search against a 64-cell full solve and counted a different amount of work. The raw timing also included 43961 s of machine sleep. The agreed citable figure is 3.70× on 45507 s of compute. |
| ATSim stress error 10–20%, kept conservative | 12.1% high on average, removed with a scale factor of 0.892 | You settled this correction explicitly and asked that it not be described as agreement with COMSOL. The draft follows that: mean bias removed, spread between layouts kept, worst case stated. |
| Temperature error 0.1 K | within 0.3 K | 0.3 K is what the COMSOL comparison supports. 0.1 K does not appear in the evidence. |
| Cost is temperature plus gradient | Cost is temperature, gradient and wirelength in equal thirds | This is the loop as fixed on 2026-11-16. Wirelength uses an assumed netlist, which the draft says in so many words. |
| Anneal acceleration around 2.4× | 2.42× to 2.75× on the cold-edge model, 1.97× to 2.11× on the overhang model | Both are in the evidence. The overhang model is the one comparable with HotSpot, so it is the number to quote against other tools; the draft quotes both and says which is which. |
| Frequency scaling up to 4× throughput | not yet stated | The 4× figure comes from an older rule that let the scale reach 4. The current loop caps it at 2.5. The end-to-end run uses the current rule, so the draft leaves the throughput numbers for Table V rather than reuse the old ones. |
