# Adding contacTrees to a BEAST xml

In this example, I create a minimal XML file for a linguistic phylogeny using BEASTling on readily available data. The XML is constructed with per-concept rate variation and otherwise minimal, to ensure that the XML file contains separate likelihoods for different concepts (which is necessary for adding contacTrees) but is otherwise relatively easy to work with.

## Preparation

The best way to use contacTrees is to add the necessary BEAST2 objects to an existing configuration file. The focus of this example is on the contacTrees, not on the other elements of configuration, so I use a standard dataset (Grollemund et al.'s Bantu data) with a standard tool (BEASTling) to generate a base XML, which I will then extend to be a contacTrees model below.

### Data

Of the four large language families with extensive phylogenetic studies, Bantu stands out for haviing an easily available data set (Indo-European's picture is much more complicated) which does not include forms from other language families (like the big Austronesian dataset, ABVD, does) and dialect chains seem to be a lesser factor in its formation (as opposed to Sino-Tibetan).

I therefore base this exapmle on the Bantu dataset from

> Grollemund, Rebecca & Branford, Simon & Bostoen, Koen & Meade, Andrew & Venditti, Chris & Pagel, Mark. 2015. Bantu expansion shows that habitat alters the route and pace of human dispersals. Proceedings of the National Academy of Sciences 112(43). 13296–13301. (doi:10.1073/pnas.1503793112)

in the CLDF version available from https://github.com/lexibank/grollemundbantu ([v1.0rc6](https://github.com/lexibank/grollemundbantu/releases/tag/v1.0rc6)). The data set is available under the CC-BY-NC 4.0 International license, so I have included a copy of it in [this folder](cldf/).

### BEASTling

To create the XML file, I use the python helper [BEASTling](https://github.com/lmaurits/BEASTling/), which from a small description file can generate a full XML file when given a CLDF dataset.

> Maurits, Luke & Forkel, Robert & Kaiping, Gereon A. & Atkinson, Quentin D. 2017. BEASTling: a software tool for linguistic phylogenetics using BEAST 2. PLOS ONE 12(8). e0180908. (doi:10.1371/journal.pone.0180908)

I use the most recent published version of BEASTling, v1.5.1 (installed using `pip install beastling`) with the [following configuration file](bantu-beastling.conf).

> ```
> [admin]
> basename=bantu-beastling
> [model bantu-vocabulary]
> model = binaryctmc
> data = cldf/cldf-metadata.json
> rate_variation=True
> ```
> — bantu-beastling.conf

Generating the BEAST2 XML file using `beastling bantu-beastling.conf` creates [the basic XML file in this folder](bantu-beastling.xml), which will serve as the basis for adding contacTrees in the following steps.

(Due to [a bug in BEASTling](https://github.com/lmaurits/BEASTling/issues/261), I had to substitute all the `,_` in the resulting XML file with just `_` before the file was accepted by BEAST2.)

## Adding contacTrees

### ConversionGraph instead of Tree

Where a standard phylogenetic analysis works with a single tree, contacTrees assumes a network composed of a backbone tree with conversion edges. The first adjustment we need to make is therefore to replace the tree stateNode with a `contactrees.ConversionGraph`, by just replacing or adding the `spec` attribute.

```
-      <tree id="Tree.t:beastlingTree" name="stateNode">
+      <tree id="Tree.t:beastlingTree" name="stateNode" spec="contactrees.ConversionGraph">
```

### Model parameters

The contacTrees model has a few additional parameters compared to a binary tree phylogenetic inference. These need to also be added as state nodes.

The first of these parameters is the conversion rate, which governs the number of conversion edges (i.e. discrete language contact events) in the tree. Here I use the parameterization where instead of a conversion rate, the expected number of conversion edges in the whole tree is specified. To do this, contacTrees needs to adjust the conversionRate to give the expected number of edges, so it needs to know the exact shape of the network (or at least the backbone tree).

    <parameter id="expectedConversions" name="stateNode" estimate="false" value="0.25"/>
    <stateNode id="conversionRate" spec="contactrees.model.ConversionRate" expectedConversions="@expectedConversions" linearContactGrowth="true" network="@Tree.t:beastlingTree"/>

Each concept forms a “block” that is either completely inherited along a conversion edge or not at all. The state needs to track for each block which path it takes, so each Block is practically a vector of binary borrow/inherit choices, with logic attached to it that deals with the adding and removing of edges. I construct the blocks based on a `plate` of concepts (which is generated elsewhere by BEASTling because of the rate variation, I only need to copy it), which I will re-use below. (I could avoid the repeated lists of concepts using some XML entity magic, but that would distract from the content of this tutorial.)

      <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
        <stateNode spec="contactrees.Block" id="$(concept)"/>
      </plate>

The “borrowing probability” (or, more generally, the probability `pMove` of a block moving over a conversion edge, instead of being inherited) is another model parameter, which we infer to be between 0 and 0.4 – the theoretical maximum is 0.5, because the only thing that internally distinguishes borrowing from inheritance is that borrowing is the minority path, while inheritance is the majority path.

    <parameter id="pMove" name="stateNode" estimate="true" value="0.1" lower="0.0" upper="0.4"/>

### Priors

The tree prior (in my case a simple YuleModel) needs to we wrapped inside an ACGDistribution. So instead of

        <distribution id="YuleModel.t:beastlingTree" tree="@Tree.t:beastlingTree" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@birthRate.t:beastlingTree" />

I have

        <distribution id="ACGTreePrior" spec="contactrees.model.ACGDistribution" network="@Tree.t:beastlingTree" conversionRate="@conversionRate" linearContactGrowth="true" upperConvCountBound="50">
          <distribution name="cfModel" id="YuleModel.t:beastlingTree" tree="@Tree.t:beastlingTree" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@birthRate.t:beastlingTree" />
        </distribution>


The each “block” may or may not move along each conversion edge. In total, the proportion of blocks that move should follow the probability `pMove`, so the following prior implements essentially a binomial distribution over all blocks and conversion edges.

      <distribution id="ConvMovePrior" spec="contactrees.model.ConversionMovePrior" network="@Tree.t:beastlingTree" pMove="@pMove">
        <blockSet spec="contactrees.BlockSet" id="allBlocks" network="@Tree.t:beastlingTree">
          <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
            <block idref="$(concept)"/>
          </plate>
        </blockSet>
      </distribution>

### Likelihoods

Each individual concept likelihood remains a tree likelihood. Hovever, the the separate word trees, or more precisely, the concept trees may differ and need to be extracted from the `contactrees.ConversionGraph` using `contactrees.MarginalTree`s. So each likelihood, instead of the attribute

    tree="@Tree.t:beastlingTree"
    
needs to contain the corresponding marginal tree for the specific block

    <tree id="marginalTree.CONCEPT" spec="contactrees.MarginalTree" network="@Tree.t:beastlingTree" block="@CONCEPT" branchRateModel="@StrictClockModel.c:default" />
    
(This can be solved with regular-expression-capable search-and-replace and is a bit cumbersome otherwise.)

### Operators

The operators for adding, removing, and changing conversion edges need to be added.

    <operator id="AddRemoveConversion.t" spec="contactrees.operators.AddRemoveConversionGibbs" weight="50.0" acg="@Tree.t:beastlingTree" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGTreePrior">
      <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
        <treeLikelihood idref="treeLikelihood.$(concept)"/>
      </plate>
    </operator>
    
    <operator id="GibbsSampleMovesPerConversion.t" spec="contactrees.operators.GibbsSampleMovesPerConversion" weight="10.0" acg="@Tree.t:beastlingTree" pMove="@pMove" blockSet="@allBlocks" mcmcmc="true">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>

    <operator id="ConvertedEdgeSlide.t" spec="contactrees.operators.ConvertedEdgeSlide" acg="@Tree.t:beastlingTree" apertureSize="0.2" weight="15.0"/>
    <operator id="ConvertedEdgeFlip.t" spec="contactrees.operators.ConvertedEdgeFlip" acg="@Tree.t:beastlingTree" weight="1.0"/>
    <operator id="ConversionSplit.t" spec="contactrees.operators.ConversionSplit" acg="@Tree.t:beastlingTree" weight="1.0"
            blockSet="@allBlocks" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" flip="false" pMove="@pMove"/>
    <operator id="ConvertedEdgeHop.source" spec="contactrees.operators.ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" sourceOnly="true" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="2.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.source.narrow" spec="contactrees.operators.ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" sourceOnly="true" nClosestRelatives="3" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="6.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.narrow" spec="contactrees.operators.ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" blockSet="@allBlocks" nClosestRelatives="3" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="2.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>

And the existing tree operators (everything that changes the language tree) need to be replaced by the corresponding contacTrees operators. These adapted tree operators ensure that conversion edges do not become invalid when the language tree is changed.

```
-    <operator id="SubtreeSlide.t:beastlingTree" spec="SubtreeSlide" tree="@Tree.t:beastlingTree" markclades="true" weight="15.0" />
-    <operator id="narrow.t:beastlingTree" spec="Exchange" tree="@Tree.t:beastlingTree" markclades="true" weight="15.0" />
-    <operator id="wide.t:beastlingTree" isNarrow="false" spec="Exchange" tree="@Tree.t:beastlingTree" markclades="true" weight="3.0" />
-    <operator id="WilsonBalding.t:beastlingTree" spec="WilsonBalding" tree="@Tree.t:beastlingTree" markclades="true" weight="3.0" />
-    <operator id="UniformOperator.t:beastlingTree" spec="Uniform" tree="@Tree.t:beastlingTree" weight="30.0" />
-    <operator id="treeScaler.t:beastlingTree" scaleFactor="0.5" spec="ScaleOperator" tree="@Tree.t:beastlingTree" weight="3.0" />
-    <operator id="treeRootScaler.t:beastlingTree" scaleFactor="0.5" spec="ScaleOperator" tree="@Tree.t:beastlingTree" rootOnly="true" weight="3.0" />
-    <operator id="UpDown" spec="UpDownOperator" scaleFactor="0.5" weight="3.0">
-      <tree idref="Tree.t:beastlingTree" name="up" />
-      <parameter idref="birthRate.t:beastlingTree" name="down" />
-    </operator>
+    <operator id="CFWilsonBalding" spec="contactrees.operators.CFWilsonBalding" acg="@Tree.t:beastlingTree" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" includeRoot="false" alpha="0.1" networkPrior="@ACGTreePrior" weight="1.0"/>
+    <operator id="CFNarrowExchange" spec="contactrees.operators.CFSubtreeExchange" acg="@Tree.t:beastlingTree" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="true" networkPrior="@ACGTreePrior" weight="10.0"/>
+    <operator id="CFWideExchange" spec="contactrees.operators.CFSubtreeExchange" acg="@Tree.t:beastlingTree" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="false" networkPrior="@ACGTreePrior" weight="1.0"/>
+    <operator id="ACGscaler" spec="contactrees.operators.ACGScaler" acg="@Tree.t:beastlingTree" scaleFactor="0.9" weight="10.0" parameterInverse="@clockRate.c:default"/>
+    <operator id="ACGscaler.rootOnly" spec="contactrees.operators.ACGScaler" acg="@Tree.t:beastlingTree" scaleFactor="0.75" weight="1.0" rootOnly="true"/>
+    <operator id="CFUniform" spec="contactrees.operators.CFUniform" acg="@Tree.t:beastlingTree" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" scaleFactor="0.9" networkPrior="@ACGTreePrior" weight="28.0"/>
```

### Logging
The logging of contacTrees results is very similar to a standard phylogenetic analysis, but in order to log the conversion edges we need to use the ACGWithMetaDataLogger. Similar to a tree-logger, this writes newick strings into a nexus file, but it uses the extended Newick format to include conversion edges.

```
-      <log id="TreeLoggerWithMetaData" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:beastlingTree" dp="4" />
+      <log id="ACGLoggerWithMetaData" spec="contactrees.ACGWithMetaDataLogger" network="@Tree.t:beastlingTree" blockSet="@allBlocks" branchratemodel="@StrictClockModel.c:default"/>
```
