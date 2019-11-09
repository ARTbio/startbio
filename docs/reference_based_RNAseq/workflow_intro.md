# Galaxy Workflows

At this point, you should be more familiar with

- importing and manipulating datasets in Galaxy
- using tools in single consecutive steps
- visualising the metadata associated to these steps as well as the results.


However, this is only the tip of the Galaxy.

Indeed, as you may have noticed, histories can become very complicated with a lot of
datasets whose origin and purpose is not so easy to remember after a while (shorter that
you may believe).

Actually, the best way to preserve an analysis is to get it completely scripted in a
computational workflow.

This is where you find the Galaxy workflows !

Galaxy workflow can be extracted from an history or built from scratch using the
Galaxy workflow editor (Menu `worflows`).

A workflow can be replayed at any time to regenerate an analysis. Importantly, they can be
exported as a `.ga` file and imported in another Galaxy server. Provided that this new
server has the input data and the tools specified by the workflow, the exact same analysis
will be generated.

Take home message: "advanced Galaxy users use workflows, to capture their work and make
convincing, transparent and re-usable their computational protocols"

In the next and last section, you will test 2 workflows that are available in your
Galaxy server and recapitulate most of the analyses you have performed today.