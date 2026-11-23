# Sections III and IV — Architecture and Evaluation

*DAC regular paper, the two sections this author owns. Target: about 3.5 pages of body together, III roughly 1.6, IV roughly 1.9. Figures are named where they go; the captions say which evidence file they come from. Numbers are taken from `论文证据/填充场叠加/README.md` and from the twelve `report.json` / `trace.json` files under `论文证据/闭环全流程_只收更优/runs/`, read on 2026-11-17. Where this text disagrees with `DAC_draft.md`, the disagreement is listed at the end and the evidence wins.*

---

## III. The Co-Design Loop

The loop has two halves, and the split between them is the design decision the rest of this section exists to justify. The inner loop moves chiplets. It proposes a placement, prices it, and accepts or rejects it, and it does this many times, so the price of one step has to be cheap and it has to stay correct when the material between the chiplets is not silicon. The outer loop does not move anything. Once a floorplan is frozen it checks the temperature and the stress against the limits set for the part, and then it spends whatever margin is left by raising frequency. Fig. 1 draws the two halves and marks the boundary.

*[Fig. 1: one wide figure. Left, the front end: a workload and a catalog of CPU, GPU, CIM and HBM parts, four power tools, HBM power set from the CPU and GPU count. Middle, the inner loop: propose one of five moves, price it — by superposition when the fill is silicon, by eight correction steps on the starting factorization otherwise — accept or reject on the cost J. Right of a vertical dashed line, the outer loop: the frozen floorplan goes to the stress solver once, then per-chiplet frequency is raised until the temperature or stress gate stops it. Label the dashed line "floorplan frozen; nothing to the right of it moves a chiplet."]*

### III-A. What the loop takes in

The input is a list of chiplets and a task for each. CPU power comes from gem5 traces through McPAT, GPU power from Accel-Sim through AccelWattch, CIM power from AccelForge, and the HBM power from a channel formula that grows with the number of CPUs and GPUs, so the memory power follows the compute rather than being typed in by hand. Parts whose datasheet gives no die area take the area of the closest commercial part of the same class, and the catalog records which one was used.

The starting floorplan comes from three constructors — a B*-tree packer, a Monte-Carlo tree search and an L-BFGS relaxation — run at the same seed, snapped onto the 1.125 mm grid and pushed apart until every pair has at least one cell of gap. The three are scored and the lowest wins. We do not start from a hand-placed shelf, because the comparison in Section IV is meant to be against a start the loop can produce on its own.

### III-B. The cost, and what it deliberately leaves out

The inner loop minimizes

```
J = ( (T_max - T_amb) / dT_0  +  g / g_0  +  L / L_0 ) / 3
```

where `g` is the largest temperature gradient on the chip layer, `L` is the wirelength per wire, and the three divisors are the values at the starting floorplan, so the three terms start equal and none of them can dominate by units alone. Lower is better. Wirelength uses the center-to-center Manhattan distance on a connection matrix that is fixed before the anneal: 1024 each way between every CPU or GPU and the HBM, 128 each way between every CIM and the first CPU. The matrix is assumed, not routed. Freezing it is what makes the initial and the final floorplan comparable; a router rerun between them would move `L` for reasons that have nothing to do with placement.

Two things are absent from `J` on purpose. Stress is absent because the stress field is not the sum of per-chiplet fields: the constraint that produces it depends on the whole temperature distribution at once, so it cannot be cached and reused the way a temperature field can, and one stress solve costs about as much as thirty thermal solves. Temperature and stress limits are also absent. A candidate is not rejected for being hot. It is rejected for being worse than the current floorplan on `J`. The limits are applied once, after the floorplan freezes, which is Section III-D.

### III-C. Moves and acceptance

A step does one of five things. Swapping a CPU or GPU for another part from the catalog is drawn 10% of the time. CIM chiplets are never swapped, and the HBM power is recomputed whenever the CPU or GPU count changes, which keeps the front end and the search consistent with each other. The other 90% is split four ways: rotation takes a fixed tenth of it, a small nudge takes a share that grows linearly from 0.10 to 0.70 as the anneal proceeds, and the rest is divided equally between a long jump and a swap of two positions. The search therefore moves coarsely at the start and finely at the end.

A move that lowers `J` or leaves it unchanged is always taken. In the runs reported here a move that raises `J` is never taken: we dropped the probabilistic acceptance of worse floorplans after it lengthened the search without improving the floorplans it found. The best floorplan seen is what is kept, and it is rescored exactly before it is handed to the outer loop.

### III-D. Pricing a move when the material map moves with it

This is the part of the architecture that exists only because the fill is not silicon, and it is the reason the method has the shape it does.

For a linear conduction problem the temperature field of several sources is the sum of the fields of each source alone. If the material map never changes, the field of chiplet `i` depends only on its own power and its own position, so we compute it once, store it, and rebuild the package field by adding stored fields. A move touches one or two chiplets, only those fields are recomputed, and the rest are read back. That is superposition, and it is exact for as long as the conductance matrix is fixed. It is the right tool when the gaps are filled with silicon, because then every cell has the same conductivity and moving a chiplet changes nothing in the matrix.

It stops being exact the moment the fill differs from silicon. Moving a chiplet moves the silicon cells and the fill cells with it, so the matrix changes. Pricing the move on the matrix from the starting layout means pricing it on a material map the package no longer has. The error is not something that accumulates politely over a long search. On a four-chiplet floorplan with epoxy in the gaps, rotating one chiplet at step 1 put the quoted peak 18.6 K above the true peak and one cell 29.2 K above it, and by step 10 the worst cell was 31.4 K off (Table I). For a search that accepts on a tenth of a kelvin, a price that is wrong by tens of kelvin at the first move is not a price.

*[Table I: ten steps on one four-chiplet floorplan, epoxy fill, both methods walked over identical proposed layouts. Columns: step, move, quoted peak with the material map frozen at the start, true peak with the map rebuilt, peak error, worst cell. Source: 论文证据/填充场叠加/README.md, the ten-step table. Show steps 1, 3, 6 and 10 in the paper; the full ten rows go in the repository.]*

| Step | Move | Quoted °C | True °C | Peak error K | Worst cell K |
| --- | --- | --- | --- | --- | --- |
| 1 | rotate GPU | 88.29 | 69.71 | +18.59 | 29.18 |
| 3 | swap HBM3 with Penryn | 88.29 | 64.53 | +23.76 | 29.18 |
| 6 | swap, rejected | 96.75 | 69.99 | +26.76 | 29.15 |
| 10 | jump GPU | 90.17 | 64.22 | +25.95 | 31.39 |

The errors sit on the cells that were silicon at the start and are fill after the move, which is exactly where a stale material map must be wrong. Rebuilding and refactoring the matrix at every step removes the error and removes the speedup with it: refactoring is where the time goes, by more than an order of magnitude over filling the matrix or back-substituting through it.

So we keep the factorization and correct for the matrix we no longer have. Let `A0` be the matrix at the start, already factored, and `A1` the matrix after the move. Back-substituting the new power through the old factor gives a trial field `T*`. The power that would actually produce `T*` on the new map is `P* = A1 T*`, which is one multiply, and the gap `P − P*` is the power the old factor failed to account for. Feeding that gap back through the same factor and repeating is a fixed-point iteration whose fixed point is the exact field of `A1`. We run it eight times and stop. There is no residual test and no refactor when a move is accepted, so the cost of a step does not depend on what the search decides.

Eight is a fixed count for two reasons. The residual after eight steps is far below the threshold the anneal accepts on, so more iterations would buy accuracy the search cannot use. And a fixed count keeps the cost of a step predictable, which matters once the same loop is run over packages from 5 to 11 chiplets and the two times have to be compared. There is no preconditioner and no rank update; the iteration is the whole of the method.

The correction needs the changed matrix to be a perturbation of one that is already a reasonable approximation. A fill as insulating as air breaks that. With air in the gaps, moving a GPU onto cells the starting map still calls air put the quoted peak 356 K above the true one, because the trial field starts from a cell that cannot conduct. We state the limit plainly: the correction is for mold compounds and underfills, not for air, and we do not run it on air.

### III-E. The stack the price is computed on

The pricing model is the HotSpot stack: a 150 µm chip layer, 20 µm of interface material at 4 W/m·K, a 1 mm copper spreader and a 6.9 mm copper sink, both 96 mm on a side over a 72 mm chip frame, convection at the sink's value, ambient at 45 °C. Chip cells are silicon at 130 W/m·K. Gap cells take the fill conductivity, which is an input — silicon, epoxy at 0.8 W/m·K, underfill at 1.0, or a named value. The 12 mm of spreader that overhangs the chip frame is modeled as eight extra nodes rather than cut off, because cutting it off raises the whole board by about 3 K against HotSpot and would make every comparison with HotSpot meaningless. The stress solver uses the same stack plus a 1 mm substrate and a 100 µm interposer, and the same convection.

### III-F. The outer loop

Once the floorplan is frozen it is not moved again. The stress solver runs on the 2 mm grid and reports the hottest die cell and the highest von Mises stress inside a normal chiplet; HBM is not counted as a normal chiplet for the stress gate. Frequency is then adjusted per CPU and per GPU, in steps of 0.1, scaling dynamic power only. The search raises scales for as long as a guide estimate stays under 90 °C, checks the true result, and walks the worst offender back down if the gate is missed. Stress above 200 MPa but within 220 MPa gets a short repair of at most six tries; stress above 220 MPa is left as it is and marked as failing the limit, so a failing floorplan stays in the results instead of disappearing.

The split is the point. The inner loop answers where the chiplets should go, and it answers that under whatever fill the package actually uses. The outer loop answers, for a floorplan that is already fixed, how far each chiplet can be pushed before the package crosses a limit. The limit is an input. The same loop serves a part with a tight limit and a part with a loose one, and the frequency scaling lands at a different point for each.

---

## IV. Evaluation

### IV-A. Setup

All anneals use seed 22 for the method comparisons and seed 11 for the end-to-end runs, and both seeds are fixed before looking at the results. The method comparisons in IV-B and IV-C price on the Python thermal field, not on HotSpot and not on the stress solver, so their absolute temperatures are not comparable with the end-to-end numbers in IV-D; each table is self-contained and says which field it uses. The end-to-end runs report three temperatures for the same floorplan — the Python field used for pricing, HotSpot, and the stress solver on its 2 mm grid — precisely so that the gap between them is visible rather than hidden behind one column.

The packages are 5, 6, 7, 8, 9 and 11 chiplets, each run with silicon fill and with epoxy fill, twelve runs in all. Every run starts from the same constructor output, takes 40 annealing steps, and accepts only moves that do not raise `J`.

### IV-B. What the fill changes

Before asking how fast the search is, it is worth seeing whether the fill changes the answer at all. If it does not, the correction in Section III-D is unnecessary.

*[Table II: share of package power leaving by each path, five fills, both secondary exits open. Source: 论文证据/填充场叠加/README.md, the percentage table; the underlying powers sum to 154.79 W within 1e-4 W. Fig. 2 draws the same rows as stacked bars.]*

| Fill | k (W/m·K) | Spreader % | Underside % | Edge % | Secondary paths % |
| --- | --- | --- | --- | --- | --- |
| Air | 0.028 | 89.62 | 7.33 | 3.05 | 10.38 |
| Organic | 0.5 | 77.11 | 10.38 | 12.50 | 22.89 |
| Epoxy | 0.8 | 73.56 | 10.34 | 16.10 | 26.44 |
| Underfill | 1.0 | 71.65 | 10.24 | 18.11 | 28.35 |
| Silicon | 130 | 30.96 | 5.04 | 64.00 | 69.04 |

*[Fig. 2: stacked bars of heat flow through spreader, underside and edge, for the five fills. Source: 论文证据/填充场叠加/图/散热路径占比百分比.png.]*

With air in the gaps, nine tenths of the power leaves through the spreader. With silicon in the gaps, the spreader carries under a third and the edges carry 64%. That is 59 percentage points of the power moved off the main path by a choice that most floorplanning tools do not expose, because they fill the gaps with the chip-layer conductivity and call it done. A search that assumes silicon fill is optimizing a package whose heat leaves by a different route from the one being built, and it will prefer different floorplans for it. The fill has to be an input.

The same table is why air is excluded from the correction rather than merely discouraged. The secondary paths carry 10% of the power with air and 69% with silicon; the matrix for an air-filled package is not a small correction of the silicon matrix, and Section III-D already showed what the iteration does when it is started that far from the solution.

### IV-C. The correction against refactoring every step

The comparison that matters is not between two anneals, because the two accept on their own quoted peaks and diverge onto different floorplans. It is between the price the iteration quoted for the floorplan it chose and an exact refactor of that same floorplan. That is the error a search would act on.

*[Table III: six packages, epoxy fill, cold package edge, 200 steps, seed 22. Time is the anneal only; the exact rescore of the chosen floorplan is timed separately and excluded. Peak error is the quoted peak minus the rescored peak, so a negative value means the search was quoted a temperature below the truth. Source: 论文证据/填充场叠加/README.md, the 5-to-11-chiplet table.]*

| Chiplets | Power W | Refactor s | Corrected s | Speedup | Worst cell K | Peak error K |
| --- | --- | --- | --- | --- | --- | --- |
| 5 | 68.89 | 119.0 | 43.8 | 2.72× | 1.003 | −0.114 |
| 6 | 140.69 | 134.8 | 55.6 | 2.42× | 2.404 | −0.777 |
| 7 | 71.08 | 153.8 | 55.9 | 2.75× | 1.094 | −0.180 |
| 8 | 138.81 | 143.2 | 54.3 | 2.64× | 2.180 | +0.081 |
| 9 | 138.84 | 136.9 | 53.4 | 2.56× | 1.192 | −0.637 |
| 11 | 209.96 | 136.1 | 51.2 | 2.66× | 2.534 | −0.203 |

The speedup sits between 2.42× and 2.75× and does not grow with the chiplet count, which is what the cost model predicts: refactoring dominates the time per step, both sides pay for the matrix fill and the eight multiplies, and the factorization is paid once. The peak error stays within 0.8 K, and it is negative in five of the six cases, so the search is quoted a slightly optimistic temperature and occasionally accepts a floorplan it would have rejected on the true value. That is the behavior to design for rather than to average away, and it is why the chosen floorplan is rescored exactly before it leaves the inner loop. The exact rescore is one solve. It is part of the method, not an optional check.

*[Fig. 3: paired bars of anneal time for the six packages, refactor against eight correction steps, with the peak error printed over each corrected bar. Source: 论文证据/填充场叠加/图/加速对比_用时.png and 图/加速对比_误差.png.]*

Repeating the comparison on the overhang model, which is the one comparable with HotSpot, gives a smaller speedup and a much smaller error.

*[Table IV: the same correction on the overhang model, epoxy fill, packages with matched chiplet area. Source: 论文证据/填充场叠加/README.md, the overhang table. Fig. 4 draws it.]*

| Chiplets | Power W | Refactor s | Corrected s | Speedup | Worst cell K | Peak error K |
| --- | --- | --- | --- | --- | --- | --- |
| 5 | 199.81 | 81.0 | 38.4 | 2.11× | 0.026 | 0.000 |
| 8 | 402.42 | 83.2 | 40.5 | 2.06× | 0.051 | 0.000 |
| 11 | 482.22 | 84.0 | 42.7 | 1.97× | 0.073 | 0.000 |

*[Fig. 4: time and cell error on the overhang model at 5, 8 and 11 chiplets. Source: 论文证据/填充场叠加/图/外伸环_加速对比.png.]*

The speedup falls to between 1.97× and 2.11× because the overhang matrix factors faster, so there is less to save. The cell error stays under 0.08 K and the quoted peak matches the rescore to a thousandth of a kelvin. This is the number to quote against tools that use HotSpot, and the cold-edge number is the one to quote for the five-layer model. They answer different questions and we report both.

### IV-D. The full loop

The twelve end-to-end runs use the loop as specified: 40 steps, seed 11, only non-worsening moves accepted, epoxy and silicon fill on each package size. Temperature here is the stress solver's reading on the 2 mm grid, the hottest chip by mean die temperature, because that is the reading the gate sees. The Python field used for pricing runs higher and is given in the repository alongside these numbers; the two are not averaged.

*[Table V: initial against best floorplan, then the same floorplan after frequency scaling. ΔT and Δσ are best minus initial, so negative is an improvement. The frequency column is the lowest scale applied to any CPU or GPU; every other one reached 2.5. Source: 论文证据/闭环全流程_只收更优/runs/*/report.json.]*

| Case | T init °C | T best °C | ΔT K | σ init MPa | σ best MPa | T after freq °C | σ after freq MPa | Min scale | Gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5, epoxy | 60.97 | 60.62 | −0.34 | 47.1 | 50.7 | 75.4 | 79.8 | 2.5 | pass |
| 5, silicon | 60.88 | 60.04 | −0.84 | 47.2 | 46.2 | 74.4 | 74.5 | 2.5 | pass |
| 6, epoxy | 67.29 | 65.09 | −2.20 | 85.0 | 75.0 | 80.8 | 111.4 | 2.5 | pass |
| 6, silicon | 64.78 | 64.57 | −0.20 | 85.8 | 77.7 | 78.3 | 114.8 | 2.5 | pass |
| 7, epoxy | 65.71 | 64.39 | −1.32 | 79.1 | 68.7 | 80.8 | 108.9 | 2.5 | pass |
| 7, silicon | 65.69 | 63.00 | −2.69 | 65.9 | 60.1 | 80.1 | 103.2 | 2.5 | pass |
| 8, epoxy | 72.42 | 72.83 | +0.41 | 96.6 | 91.7 | 88.7 | 153.8 | 2.5 | pass |
| 8, silicon | 72.17 | 72.57 | +0.40 | 88.6 | 84.6 | 88.7 | 140.5 | 2.5 | pass |
| 9, epoxy | 70.11 | 70.18 | +0.07 | 104.8 | 94.0 | 88.4 | 147.5 | 2.5 | pass |
| 9, silicon | 70.00 | 70.07 | +0.07 | 112.3 | 96.5 | 88.4 | 152.1 | 2.5 | pass |
| 11, epoxy | 74.21 | 73.51 | −0.70 | 105.7 | 102.5 | 89.2 | 152.5 | 2.2 | pass |
| 11, silicon | 73.98 | 73.21 | −0.78 | 112.0 | 108.5 | 89.1 | 161.6 | 2.2 | pass |

*[Fig. 5: for each of the six package sizes, three bars — initial, best, after frequency scaling — of the hottest-chip temperature, epoxy and silicon side by side. Source: the same report files.]*

Three things come out of the table, and none of them requires the best case to make the point.

Placement and frequency do different jobs, and the table separates them. The largest temperature drop from placement is 2.7 K, on the 7-chiplet silicon package. Frequency scaling then spends the margin that placement created: every package ends between 74 °C and 89 °C, up against the 90 °C gate, with the stress still under 200 MPa in all twelve. The gain from frequency is not a lower temperature. It is the scale, and on ten of the twelve packages every CPU and GPU reached 2.5, the top of the range we searched. Only the 11-chiplet packages ran out of margin first, and only on one chiplet, which stopped at 2.2 while the others reached 2.5. Reporting a single package-wide scale would have hidden that.

The fill changes the stress more than it changes the temperature. At 11 chiplets the epoxy package freezes at 102.5 MPa and the silicon package at 108.5 MPa, from starting points of 105.7 and 112.0. The temperature gap between the two fills is under 0.4 K on the same package. A cost function that sees only temperature would call these two packages equivalent. The stress gate does not, which is the practical argument for computing the fill rather than assuming it.

The 8-chiplet packages are the honest counterexample, and we keep them. Their cost fell by half a percent, but the temperature rose by 0.4 K, because the accepted moves traded temperature for wirelength: the wirelength per wire dropped from 35.0 mm to 34.1 mm while the gradient barely moved. An objective with three equal terms will do this whenever the start is already good on temperature, and the 8-chiplet start was. A reader who looks only at the temperature column would say the search failed. It did what it was asked. We would rather show that than drop the rows.

### IV-E. What we would not claim

The correction is not reliable for a fill as insulating as air, and IV-B shows why the gap is too wide to be called a corner case. The wirelength uses an assumed netlist rather than a routed one; it is frozen so that the comparison is fair, but it is not a routing result. The frequency search reported here stops at a scale of 2.5, and on most of these packages it stopped because it hit that end of the range rather than because it hit the thermal gate, so the throughput gain is limited by our search range and not by the package. We do not quote a throughput multiple from these runs.

---

## Notes on disagreements with DAC_draft.md

Kept here so the two files can be reconciled. The evidence wins in each case.

| Draft | This file | Why |
| --- | --- | --- |
| Drift of the correction "0.37 K to 0.85 K" | not used | The twelve end-to-end traces give a rescored drift of 0.34 K to 0.80 K on epoxy and exactly 0 on silicon, which is a different and narrower claim than the draft's. Section III-C points at the exact rescore instead of quoting a range the draft could not source. |
| "its danger is the candidate that lands on the threshold, where a 0.02 K pricing error flips the decision" | not repeated here | That observation belongs to the superposition-against-TAP comparison, which is outside these two sections. |
| Frequency range stated as 0.5 to 2.5 | stated, and flagged as binding | Table V shows ten of twelve packages pinned at 2.5, so the cap is doing the work, not the gate. IV-E says so rather than reporting the scales as thermal results. |
| Anneal "priced by superposition" throughout | superposition only when the fill is silicon | With epoxy the matrix changes and superposition is the method Section III-D shows failing by 18.6 K at step 1. The inner loop uses the eight-step correction there. |
