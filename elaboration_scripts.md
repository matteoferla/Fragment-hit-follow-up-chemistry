# Elaboration pipelines
For full details see [![Read the Docs](https://img.shields.io/readthedocs/fragment-hit-follow-up-chemistry)](https://fragment-hit-follow-up-chemistry.readthedocs.io/en/latest/fragment_elaboration_scripts.html)


## Fragmenstein

In [fragment_elaboration_scripts.fragmenstein_merge_sw_place](fragment_elaboration_scripts/fragmenstein_merge_sw_place.py)
is the old script to run Fragmenstein on a PDB file and SDF hits, merge them,
find analogues via SmallWorld and place them in the PDB.
It is now fully incorporated into Fragmenstein.

> For more [:octocat::Fragmenstein](https://github.com/matteoferla/Fragmenstein)

## Fragment network

> :construction: Clean-up

## Arthorian quest

> See [:octocat::Arthorian quest](https://github.com/matteoferla/arthorian-quest)

> :construction: Repo under construction

## Filtering

> :construction: Clean-up

## CoPriNet

> See [:octocat::CoPriNet](https://github.com/oxpig/CoPriNet)

Predict the price for a compound.

> :construction: hunt down API code


## OE ROCS

This task has three parts:

1. Subsetting of Enamine BB, see below.
2. Conformer generation with OE Omega.
3. FastROCS with the conformers, or OE ROCS GUI derp around.

In [:octocat:: Functional-subsets-of-Enamine-BB](https://github.com/matteoferla/Functional-subsets-of-Enamine-BB)
Enamine BB is split into functional subsets.
Due to size, the OE Omega conformers are not included.

These are generated with [oe_confomer_gen.py](fragment_elaboration_scripts/oe_confomer_gen.py),
which is a bit more than a wrapper for OE Omega.

> :construction: Copy 2 & 3 off cluster

> :construction: Hunt down the fluorescent subset!
