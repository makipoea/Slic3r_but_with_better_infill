See the [project homepage](https://slic3r.org/) at slic3r.org for more information.

### Directory structure

* `package/`: the scripts used for packaging the executables
* `src/`: the C++ source of the `slic3r` executable and the CMake definition file for compiling it
* `src/GUI`: The C++ GUI.
* `src/test`: New test suite for libslic3r and the GUI. Implemented with [Catch2](https://github.com/catchorg/Catch2)
* `t/`: the test suite (deprecated)
* `utils/`: various useful scripts
* `xs/src/libslic3r/`: C++ sources for libslic3r
* `xs/t/`: test suite for libslic3r (deprecated)
* `xs/xsp/`: bindings for calling libslic3r from Perl (XS) (deprecated)
* `xs/src/libslic3r/Fill/`: here, I hope are the most modification
  -> Infill added : Nothing

Pour compiler il faut juste que tu install perl (si possible depuis le terminale) que tu execute dans le dossier source perl Build.pl // pour des compilation partiels notament aprés des modifs des fichiers xs/src/libslic3r/Fill il ne suffit pas de recompiler les fichiers objets en .obj mais on peut recompiler uniquement xs avec "perl Build -xs" ce que j'ai mis 2 jours a trouver

### Acknowledgements

The main author of Slic3r is Alessandro Ranellucci (@alranel, *Sound* in IRC, [@alranel](http://twitter.com/alranel) on Twitter), who started the project in 2011.

Joseph Lenox (@lordofhyphens, *LoH* in IRC, [@LenoxPlay](http://twitter.com/LenoxPlay) on Twitter) is the current co-maintainer.

Contributions by Henrik Brix Andersen, Vojtech Bubnik, Nicolas Dandrimont, Mark Hindess, Petr Ledvina, Y. Sapir, Mike Sheldrake, Kliment Yanev and numerous others. Original manual by Gary Hodgson. Slic3r logo designed by Corey Daniels, <a href="http://www.famfamfam.com/lab/icons/silk/">Silk Icon Set</a> designed by Mark James, stl and gcode file icons designed by Akira Yasuda.

### How can I invoke Slic3r using the command line?

The command line is documented in the relevant [manual page](https://manual.slic3r.org/advanced/command-line).
