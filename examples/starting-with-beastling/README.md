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

The contacTrees model has a few additional parameters compared to a binary tree phylogenetic inference. These need to also be added as state nodes. Here I use the parameterization where instead of a conversion rate, the expected number of conversion edges in the whole tree is specified. (It then needs to distribute the rate on the entire tree, so it needs to know about the conversion graph.)

    <parameter id="expectedConversions" name="stateNode" estimate="false" value="0.25"/>
    <stateNode id="conversionRate" spec="contactrees.model.ConversionRate" expectedConversions="@expectedConversions" linearContactGrowth="true" network="@Tree.t:beastlingTree"/>

...

    <parameter id="pMove" name="stateNode" estimate="true" value="0.1" lower="0.0" upper="0.4"/>

...

      <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
        <stateNode spec="contactrees.Block" id="$(concept)"/>
      </plate>

### Priors

The tree prior (in my case a simple YuleModel) needs to we wrapped inside an ACGDistribution. So instead of

        <distribution id="YuleModel.t:beastlingTree" tree="@Tree.t:beastlingTree" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@birthRate.t:beastlingTree" />

I have

        <distribution id="ACGTreePrior" spec="contactrees.model.ACGDistribution" network="@Tree.t:beastlingTree" conversionRate="@conversionRate" linearContactGrowth="true" upperConvCountBound="50">
          <distribution name="cfModel" id="YuleModel.t:beastlingTree" tree="@Tree.t:beastlingTree" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@birthRate.t:beastlingTree" />
        </distribution>


The model tracks for each concept (“block”) which conversion edges are used or not used. I construct this based on a `plate` of concepts (which is generated elsewhere by BEASTling because of the rate variation).

      <distribution id="ConvMovePrior" spec="ConversionMovePrior" network="@Tree.t:beastlingTree" pMove="@pMove">
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

    <tree id="marginalTree.two" spec="MarginalTree" network="@Tree.t:beastlingTree" block="@two" branchRateModel="@clock.medium"/>
    
### Operators

The operators for adding, removing, and changing conversion edges need to be added.

    <operator id="AddRemoveConversion.t" spec="AddRemoveConversionGibbs" weight="50.0" acg="@Tree.t:beastlingTree" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGTreePrior">
      <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
        <treeLikelihood idref="treeLikelihood.$(concept)"/>
      </plate>
    </operator>
    
    <operator id="GibbsSampleMovesPerConversion.t" spec="GibbsSampleMovesPerConversion" weight="10.0" acg="@Tree.t:beastlingTree" pMove="@pMove" blockSet="@allBlocks" mcmcmc="true">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>

    <operator id="ConvertedEdgeSlide.t" spec="ConvertedEdgeSlide" acg="@Tree.t:beastlingTree" apertureSize="0.2" weight="15.0"/>
    <operator id="ConvertedEdgeFlip.t" spec="ConvertedEdgeFlip" acg="@Tree.t:beastlingTree" weight="1.0"/>
    <operator id="ConversionSplit.t" spec="ConversionSplit" acg="@Tree.t:beastlingTree" weight="1.0"
            blockSet="@allBlocks" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" flip="false" pMove="@pMove"/>
    <operator id="ConvertedEdgeHop.source" spec="ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" sourceOnly="true" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="2.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.source.narrow" spec="ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" sourceOnly="true" nClosestRelatives="3" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="6.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.narrow" spec="ConvertedEdgeHopGibbs" acg="@Tree.t:beastlingTree" blockSet="@allBlocks" nClosestRelatives="3" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGTreePrior" weight="2.0">
        <plate var="concept" range="animal,arm,ashes,bark,bed,belly,big,bird,bite,blood,bone,breast,burn,child,cloud,come,count,dew,die,dog,drink,ear,eat,egg,elephant,eye,face,fall,fat_oil,feather,fingernail,fire,fire-wood,fish,five,fly,four,give,goat,ground_soil,hair,head,hear,heart,horn,house,hunger,intestine,iron,kill,knee,knife,know,leaf,leg,liver,louse,man,moon,mouth,name,navel,neck,night,nose,one,person,rain,road,root,salt,sand,see,send,shame,sing,skin,sky,sleep,smoke,snake,spear,steal,stone,sun,tail,ten,three,tongue,tooth,tree,two,urine,village,vomit,walk,war,water,wind,woman">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    

