# Licence audit (source quotations; NOT legal advice, NOT a clearance)

Date: 2026-10-11. Scope: the components a CAMB-replacement build of CosmoRec.jl would use, ship or redistribute. Nothing has been bundled, published or relicensed. **This document makes no determination of legality.** It records what the sources say and where they are silent.

## 1. Original CosmoRec (the code ported, and the data tables read at runtime)

- **Source:** `cmbcheb_test/tools/CosmoRec`, git `0867690`, remote `https://github.com/cmbant/CosmoRec.git` (a fork; the upstream author is J. Chluba).
- **There is no LICENSE/LICENCE/COPYING file** in the tree (filesystem search of the top level).
- **`README:7-8`:**
  > // Author: J. Chluba (MPA & CITA & JHU & JBCA).
  > // All rights reserved. No guarantee for the correctness of the outputs is given!
- **`README:22-28`:**
  > // LICENCE:
  > // This code and all contained libraries and subroutines can be used free of charge provided that:
  > // (i) their use will be acknowledged in publications
  > // (ii) the following papers are cited:
- Papers listed at `README:30-39`: Chluba & Thomas 2010 (MNRAS 412, 748); Chluba & Sunyaev 2006 (A&A 446, 39); Chluba, Vasil & Dursi 2010 (MNRAS 407, 599); Switzer & Hirata 2008 (PRD 77, 083006); Grin & Hirata 2010 (PRD 81, 083005); Ali-Haimoud & Hirata 2010 (PRD 82, 063521); Rubino-Martin et al. 2010 (MNRAS 403, 439).
- **`README:40-44`:** further papers "are considered for citation", plus "(iii) bugs are immediately reported to Jens@Chluba.de".
- **`README:323-327`:** the recombination database (`Rec_database`) "is a large part of the distribution of 'CosmoRec' and was computed using the codes developed in collaboration with Rajat Thomas (CITA)".
- `Development/Cosmology/Recfast++/README`: no licence or copyright statement found (grep for licen/copyright/permission/distribut/rights).

**What the text grants (as written):** use free of charge, conditional on acknowledgement and citation (and bug reporting).

**What the text does NOT mention** (genuine uncertainty):
- redistribution of the source or of the data tables (`Rec_database`, `f.corr.dat`, two-photon/Raman data, HeI data);
- modification, translation, or derivative works (a Julia port);
- sublicensing, e.g. under MIT;
- whether "All rights reserved" is meant to restrict redistribution beyond what the LICENCE paragraph grants.

**Consequence for this project (policy blocker P1, unchanged):**
- **Do not bundle** any of the native data files in CosmoRec.jl or in a bridge package.
- **Do not change any public licence statement** before written permission from the copyright holder (J. Chluba; the cmbant fork maintainer may also need asking), covering:
  - (a) redistribution of the tables;
  - (b) distribution of the translated Julia code;
  - (c) the licence under which (b) may be released.
- Until then, data is read at runtime from a user-provisioned original checkout (explicit path; see `load_runmode0_tables`).

## 2. CosmoRec.jl itself

- **`CosmoRec.jl/LICENSE:1-3`:** "MIT License / Copyright (c) 2026 CosmologicalEmulators".
- The package source states it is a port of the original (`src/CosmoRec.jl:3-6`).
- 13 of 23 `src/*.jl` files carry the header "CosmoRec is (c) J. Chluba et al.; use must be acknowledged and … cited". 10 do not: CosmoRec.jl, HIDiffusionFeedback.jl, HIPDEDefine.jl, HIPDEIntegrals.jl, HIPDEProfiles.jl, HIPDEQuadPlan.jl, HIPDESolver.jl, PrefixReuse.jl, RecombinationDiffusion.jl, RHSWorkspace.jl. Adding the notice is a source edit, deferred while the suite is frozen.
- **Uncertainty:** whether an MIT licence on a translation of "all rights reserved" source is compatible with §1. This is NOT resolved here. It is a question for the copyright holder or legal counsel, and is part of P1.

## 3. CAMB (the host code; only scratch builds are made, nothing redistributed)

- `cmbcheb_test/tools/CAMB-cosmorec/LICENCE.txt:1-5`: "Copyright 2000-2025 Antony Lewis and Anthony Challinor. Portions of this software were originally based on CMBFAST …". The remaining terms were not quoted here; read the full file before any redistribution of a modified CAMB.
- A replacement design that **links a different `libCosmoRec.a` without modifying CAMB sources** creates no CAMB source changes. Distributing a modified CAMB build is out of scope.

## 4. GSL

- `cmbcheb_test/tools/camb-cosmorec-env/conda-meta/gsl-2.8-hbf7d49c_1.json`: `"license": "GPL-3.0-or-later"`.
- The native CosmoRec and the current CAMB+CosmoRec build link GSL statically (`Makefile_main:85-86`, harness build scripts).
- A Julia shim that replaces `libCosmoRec.a` would not need GSL (`GSL_LINK=""`). That is NOT yet built or verified.
- GPL obligations for any distributed binary that links GSL are a separate question; nothing is distributed today.

## 5. Julia dependencies (local inspection only, nothing installed in this chunk)

The solver packages already in the local depot (`~/.julia/packages`: OrdinaryDiffEq*, SciMLBase, PyCall, PythonCall) are distributed under their own licences, which were not audited here. Audit these before any bundling or sysimage distribution.

## 6. Required actions before any public or bundled release

1. Written permission per §1 (P1).
2. Citation list carried into the package README and the bridge outputs.
3. Licence review of the bridge dependency set (§5).
4. Decision on GPL implications if any distributed artifact links GSL (§4).
