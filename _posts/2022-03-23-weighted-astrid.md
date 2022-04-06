---
layout: post
title:  "Weighted ASTRID"
tags: phylogenetics wip
---

> We explore how to use branch support to improve ASTRID for species tree estimation, just like weighted ASTRAL.

> (This is heavily WIP and very unfinished.)

Motivated by weighted ASTRAL (wASTRAL, from the
Mirarab lab, currently in preprint), I began thinking about a weighted variant
of ASTRID. This is obviously not as easy as it seems -- quartet based species tree reconstruction
is very different from what ASTRID does, a form of "split probability" based reconstruction. I also
doubt how much mathematical rigor I can retain going from their wASTRAL to my wASTRID. We consider the problem
of weighting by support only, not including any discussion of weighting by length.


The original estimator for ASTRID is defined to be the sample average of internode distances
across the gene trees. Recall that the internode distance is the same as the number of bipartitions
that "splits" the taxa: for taxa $$a, b$$, the number of bipartitions of form $$\{a\} \cup A | \{b\} \cup B$$
in some tree. Let these splitting bipartitions in some observed gene tree (assumed
to be estimated unless otherwise stated) be denoted as $$X(a, b)$$ for taxa $$a$$ and $$b$$.
We want a good estimator for the number of splits between these two taxa in the true gene tree, that is not simply just counting what we see
in the estimated gene tree.
Let our estimator be denoted as $$d(a, b)$$ for a single gene tree, implying that the original estimator is set to $$d(a, b) = X(a, b)$$
for some gene tree. The final outputted internode distance is $$D[a, b]$$, the average of all $$d(a, b)$$. 

We reuse wASTRAL's notations when possible. Let $$A_e$$ denote the event of $$e$$ is observed in some estimated gene tree $$G$$. 
Let $$A_e^*$$ denote the event of $$e$$ appearing in the true gene tree. Let $$s(e)$$ or $$s_e$$ denote the support of a branch.
Recall that there are multiple forms of branch support for gene trees. We do not assume some particular form of support 
as it is actually an important choice. The following identity is trivial, but motivates our exposition:

$$
X(a, b) = \sum_{e \in X(a,b)} s_e^0
$$

## Choice of estimators

As it turns out there seems to be two somewhat obvious selectors of the new estimator: a dumb one and a more sophisticated one:

 - Sum of branch support (where $$s_e = 1$$ if the split is trivial)
 - Pseudo-wASTRAL like probablistic interpretation with monte-carlo approximation (yes, that's a mouthful)

We go over each approach, and it is interesting how different types of branch support might correspond well to different types ("global"
vs. "local")

### Sum of branch support between two nodes

This is less crazy than it sounds, and it works surprisingly well for non-parametric bootstrap support for the one dataset that I tested. Consider the "normal" pipeline of estimating a gene tree with non-parametric bootstrap support:

 1. Construct a best ML tree from the gene alignment
 2. Construct a number (say 100) of bootstrap replicates
 3. Map support to each split in the best ML tree: FBP support simply maps the ratio of that split appearing in the bootstrap replicates. TBE is more sophisticated, but either way we have a "global" approach to branch support.

The catch here is that **the branch support is the crude "global" approximation of the probability of that split being correct!** [citation needed, of course, but just believe it for now]. This is notably different from local approximations of branch support, for example aBayes or SH-aLRT, which I conjecture will benefit better from a more "information-theoretical" approach a la wASTRAL's definition.

This implies that the internode distance can be set as:

$$
d(a, b) = \sum_{e \in X(a,b)} s_e
$$

because recall that we are trying to estimate $$d^*(a, b)$$. Assuming that $$X(a, b)$$, the splitting edges estimated between $$a$$ and $$b$$,
are "reasonable" to a degree, $$|X(a, b)|$$ will be a worse estimator for $$d^*(a, b)$$ compared to $$d(a, b)$$. This might even have
a better property when we are using TBE. The implementation will be simple and still fast, and the problem for setting up experiments
is simply that bootstrapping is extremely costly. While things like IQTree's UFBoot2 or RAxML's RBS might help (they are the only
global branch support that I am aware of besides standard BS), it is not clear if they should be as trusted as SBS.

Now, it is interesting why intuitively I think a global branch support will work well for this simplistic approach -- one conjecture is that
local measures of branch support can work well for quartet based methods because they are literally exploring alternative
topologies between $$ab|cd$$ vs. $$ac|bd$$, etc. surrounding the branch evaluated -- in some sense quartet based already
(this is very vague, but I will hopefully make this rigorous someday), but normal branch support is decidedly split-based.

Or even there is a chance that "normal" bootstrap support is better at estimating the true internode distance compared to its intended
purpose for interpreting correctness of clades -- I don't know, and I certainly hope so, but it seems overly optimistic at this point.
Another property that we should keep in mind of is that this estimator can be computed fast: it is just the internode distance with
unequal distances between nodes.

### Somewhat information theoretic interpretation

It is useful to have an alternative to the above approach (sum of branch support) because 1) we still want to be able to use other
measures of branch support. For example aBayes looks quite good both on its own and also quite good when paired with wASTRAL, but
the current results indicate that, simply summing the branch support in this case (aBayes) seems to be quite terrible. In this case,
we want to emulate somewhat the interpretation brought by wASTRAL, namely we want something in this form:

$$
P(A_e^*|A_e) = s(e) + (1 - s(e)) B(e)
$$

where $$B(e)$$ is some background probability. To borrow wASTRAL's approach, we want the branch support to range from 0, meaning
no information at all, with some background probability being correct, but if it is 1, then it has to be correct.

The current best way that I thought of doing this is to approximate the probabilities based on the following assumption:

 - Given a set of pairwise compatible splits $$E$$ with $$s(e) = 1$$ in some estimated gene tree $$G$$, with rest of the edges, $$Bip(G) - E$$
 having $$s(e) = 0$$, then any $$S \subseteq Bip(G) - E$$ s.t. $$S \cup E$$ forms a tree is equally probable. In human language,
 an estimated gene tree that has only binary support has all its $$1$$-support branch correct, and each $$0$$-support branch should be contracted. After this process
 any random resolution of this tree is equally likely.

This assumption leads to the following Monte-carlo algorithm for approximating $$d^*(a, b)$$, the true gene tree internode distance
between $$a$$ and $$b$$:

 1. Create $$N$$ replicates of the estimated gene tree $$G$$, where each branch is randomly contracted with $$1 - s(e)$$ probability
 2. For each replicate $$G_i$$, calculate $$d_i(a, b)$$ where $$d_i$$ means the expected distance between $$a, b$$ under a random
 resolution of $$G_i$$ -- this can be somewhat easily calculated by correcting the internode distance by the polytomy degrees of each node.
 3. Return the average of all $$d_i(a, b)$$

Of course there is no proof that this a good approximation, or in any sense consistent or better than the original estimator, but
this does not seem unreasonable (hopefully). We observe that in the simple sum of bootstrap support scheme, a support of zero
will imply that we simply assume the split being incorrect, where in this case the probability of the split's existence
is dependent on the background branch support of the other branches: with low-support branches elsewhere, this branch still has a high
chance to appear "in background". This makes this algorithm less sensitive to the kind of branch support we are using -- for example
conceivably aBayes will work well as the measure of branch support here.


## Preliminary results

 > To be actually written

 The problem of course is we want gene trees that have the support we want, and if the gene tree support happens to
 be bootstrap support then reestimating the gene trees can take a really long time. Unfortunately for a split-based
 approach like ASTRID which very much relies on global information, bootstrap support might be the best but can take an awful
 lot of time to obtain. So let's wait for the bootstrap support trees to come back, so we have something beyond just S100.

<figure>
  <img src="/assets/images/wastrid/s100_prelim.png" alt="S100 Results"/>
  <figcaption><strong>Preliminary results showing nRF error rates on the S100 dataset.</strong> The non-contracted version of S100 is used.</figcaption>
</figure>

(So right now everything is promising, but certainly we need more improvements.)
