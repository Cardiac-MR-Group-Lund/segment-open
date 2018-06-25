# Contributing to segment-open

In the source code version of Segment you can incorporate your developed plugins. Please see the [Technical Manual](https://github.com/Cardiac-MR-Group-Lund/segment-open/blob/master/source/Docs/techmanual.pdf) for Segment for more details.

The developer has full rights to their plugins and may or may not release them under license terms they find appropriate as long as it does not conflict with license by Medviso AB.

We encourage users to release their plugins freely available for research purposes along with the general intentions of Segment. Freely available plugins can be incorporated in the base distribution of Segment if requested by the author and the plugin meets quality requirements.

You are not able to by yourself compile your developed plugins as part of Segment. If you are interested in a commercial distribution of your plugin please contact sales@medviso.com to discuss a partnership.

You may not distribute Segment along or bundled with your plugins, unless granted by Medviso AB.


## How to get started

 1. Start by [forking](https://help.github.com/articles/fork-a-repo/) this repository. Forking means that you create your own copy of Segment on Github.
 1. Happy hacking!
 1. When/if you want to submit your work to us, push your changes to a topic branch in your forked repository and submit a [pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) and we'll have a look at your code. If we think there is some more organization needed, we'll help.

## Things to keep in mind:

 * Found a bug? Please let us know by [opening an issue in the tracker](https://github.com/Cardiac-MR-Group-Lund/segment-open/issues)!
 * Please have a look at the coding standards in Chapter 6 of the [tech manual](../source/Docs/techmanual.pdf).
 * The [tech manual](../source/Docs/techmanual.pdf) and the Open Access [Segment paper](https://doi.org/10.1186/1471-2342-10-1) should give you an overview on how things work internally.
 * You can also look at some of the existing plugins, like [plugin_calibrate](../source/plugin_calibrate.m), [plugin_summarize](../source/plugin_summarize.m) and the basic template [plugin_template](../source/plugin_template.m) to get started.
 * Organize your code as a single file source/plugin_*pluginname*.m (to get an item in the plugin menu).
 * If you want additional files, put them in a single package directory source/+*pluginname*.
 * Use the [Matlab package mechanism](https://mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) instead of *addpath*.


# Ready to contribute?

[Take me back to the main page!](https://github.com/Cardiac-MR-Group-Lund/segment-open)
