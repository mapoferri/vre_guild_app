********
Examples
********

The best way to learn is with examples, and ``vre_template_tool`` is no exception. For this
reason, there is an example you can test and it is available at `tests/basic/`_ directory. The example is a test tool
that takes an input file, then append an username of characters in that file and finally prints the result to a second file.

.. _tests/basic/: https://github.com/inab/vre_template_tool/tree/master/tests/basic/

Default structure
~~~~~~~~~~~~~~~~~

Before delving into the details, let's first understand the directory of a test. Though it
can be modified, all VRE Tool projects have the same file structure by default, similar
to this:

::

    tests/basic/
         config.json
         in_metadata.json
         README.md
         test_VRE_RUNNER.sh


Configuration files
~~~~~~~~~~~~~~~~~~~

There are 2 configuration JSON files as inputs for the test instance. These describe the input and output files and required
arguments that need to get passed to the application. These configuration files are those that would get passed to the tool by the VRE.

The test instance will look for configuration files in:

1. ``tests/basic/config.json``
2. ``tests/basic/in_metadata.json``

config.json
===========

Defines the configurations required for by the application including parameters that need to be passed from the VRE submission form,
file and the related metadata as well as the output files that need to be produced by the application.

.. literalinclude :: ../tests/basic/config.json
   :language: json

In the arguments there 3 sets (execution, project and description) that will always be present and are provided by the VRE at the point
of submission of the to the tool. These are the name of the project that has been given in the VRE, the execution path that is
the location for where the input files are located and can be used as the working directory for the tool. Other parameters in the arguments
list are elements based on what parameters the tool requires from the user at run time.

in_metadata.json
================

List of file locations thst are used as input. The configuration names should match those that are in ``config.json`` file
defined above.

.. literalinclude :: ../tests/basic/in_metadata.json
   :language: json


Running the tool
~~~~~~~~~~~~~~~~

To run the tool it needs a config.json file and in_metadata.json file to provide de input. In this case, it uses the example
configuration files saved in ``tests/basic`` directory.

To put the tool to work, go to the project's top level directory and run:

::

    ./VRE_RUNNER --config tests/basic/config.json --in_metadata tests/basic/in_metadata.json --out_metadata out_metadata.json --log_file VRE_RUNNER.log

This command runs the tool that we've just created and you will get an output similar to this:

::

    2021-08-06 13:08:34 | INFO: 0) Unpack information from JSON
    2021-08-06 13:08:34 | INFO: 1) Instantiate and Configure Tool
    2021-08-06 13:08:34 | INFO: 2) Run Tool
    2021-08-06 13:08:34 | PROGRESS: <myApplication> execution finished successfully - FINISHED
    2021-08-06 13:08:34 | INFO: 3) Create information
    2021-08-06 13:08:34 | INFO: 4) Pack information to JSON
    2021-08-06 13:08:34 | PROGRESS: <myTool> tool successfully executed;

Results
=======

This will create a ``run000`` directory with the following contents:

::

    run000/
        goodbye.txt
        out_metadata.json
        tool.log

Now, check the files in the ``run000`` directory. You should notice that a new file has been created: ``output_metadata.json``,
with the content for the respective results, as the `run` method instructs.

.. note:: If you are wondering about the ``run000`` directory, you can see an example in the :doc:`Examples Section</examples>`.
