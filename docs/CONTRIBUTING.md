# Contributing to segment-open

 1. Start by [forking](https://help.github.com/articles/fork-a-repo/) this repository. Forking means that you create your own copy of Segment on Github.
 1. Happy hacking!
 1. When/if you want to submit your work to us, push your changes to a topic branch in your forked repository and submit a [pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) and we'll have a look at your code. If we think there is some more organization needed, we'll help.

Things to keep in mind:

 * Please read the official [Medviso AB stance on plug-ins](http://medviso.com/research/community/).
 * Found a bug? Please let us know by [opening an issue in the tracker](https://github.com/Cardiac-MR-Group-Lund/segment-open/issues)!
 * Please have a look at the coding standards in Chapter 6 of the [tech manual](source/Docs/techmanual.pdf).
 * The [tech manual](source/Docs/techmanual.pdf) and the [Segment paper](https://bmcmedimaging.biomedcentral.com/articles/10.1186/1471-2342-10-1) give you an overview on how things work internally.
 * You can also look at some of the existing plugins, like [plugin_calibrate](source/plugin_calibrate.m), [plugin_summarize](source/plugin_summarize.m) and the basic template [plugin_template](source/plugin_template.m) to get started.
 * Organize your code as a single file source/plugin_*pluginname*.m (to get an item in the plugin menu).
 * If you want additional files, put them in a single package directory source/+*pluginname*.
 * Please use the [Matlab package mechanism](https://mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) instead of *addpath*.


# Ready to contribute?

[Take me back to the main page!](https://github.com/Cardiac-MR-Group-Lund/segment-open)
