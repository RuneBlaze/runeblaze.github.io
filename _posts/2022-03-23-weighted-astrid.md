---
layout: post
title:  "Weighted ASTRID"
tags: phylogenetics wip
---

> (This is heavily WIP and very unfinished.)

> (I am showing here some personal, entirely unvetted thoughts.)

Motivated by weighted ASTRAL (wASTRAL, from the
Mirarab lab, currently in preprint and is in my opinion well written), I began thinking about a weighted variant
of ASTRID. This is obviously not as easy as it seems -- quartet based species tree reconstruction
is very different from what ASTRID does, a form of "split probability" based reconstruction. I also
doubt how much mathematical rigor I can retain going from their wASTRAL to my wASTRID.

There are two obvious choices for modifying ASTRID to wASTRID, both involving modifying
the estimator for the internode distance (defined here to be the # of edges between taxa).
The original estimator for ASTRID is defined to be the sample average of internode distances
across the gene trees. Recall that the internode distance is the same as the number of bipartitions
that "splits" the taxa: for taxa $$a, b$$, the number of bipartitions of form $$\{a\} \cup X | \{b\} \cup Y$$
in some tree.

(To be added.)

Thus this leads to the following two obvious estimators of the internode distance taking into account
the uncertainty of branches:

 - Sum of branch support between two nodes
 - Expected number of branches that "splits" the two taxa assuming that each branch $$e$$ has chance $$s_e$$
 of being correctly inferred in addition to some background probability. (The background probability is hard to describe and I really need to describe it.)

## The choice of branch support

## Preliminary results

<figure>
  <img src="/assets/images/wastrid/s100_prelim.png" alt="S100 Results"/>
  <figcaption><strong>Preliminary results showing nRF error rates on the S100 dataset.</strong> The non-contracted version of S100 is used.</figcaption>
</figure>

(So right now everything is promising, but certainly we need more improvements.)
