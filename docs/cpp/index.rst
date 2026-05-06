C++ API
=======

If Doxygen XML is available (see ``docs/Doxyfile`` and ``docs/conf.py``), the C++ API can be rendered via Breathe.

.. ifconfig:: breathe_default_project != ''

   Doxygen XML was detected, so the C++ API index is embedded below.

   .. doxygenindex::
      :project: qss

.. ifconfig:: breathe_default_project == ''

   Doxygen XML was not detected. Rebuild docs with::

      make docs

