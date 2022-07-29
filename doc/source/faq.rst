.. _faq:

Frequently Asked Questions
==========================

.. _what-version-am-i-running:

Why don't I have any absorption features in my spectrum?
--------------------------------------------------------

There are many reasons you might not have any absorption features in your
spectrum, but we'll cover a few of the basic explanations here.

 #. Your absorbers are in a different part of the spectrum than you're plotting.
    Make sure you are plotting the wavelength range where you expect to see the
    absorption by taking into account the wavelength of your absorption features
    coupled with the redshift of your dataset: :math:`\lambda_{obs} = (1 + z) \lambda_{rest}`
    To see the wavelength of specific ionic transitions, see the line list in:
    ``/trident/trident/data/line_lists/lines.txt``.

 #. Your sightlines do not have sufficient column densities of the desired
    ions to actually make an absorption feature.  Look at the total column
    density of your desired ions in your sightline by multiplying the
    density times the path length and summing it all up.  Here is an
    example for showing the total O VI column density in a ray::

        import trident
        <generate/load your ray object>
        trident.add_ion_fields(ray, ['O VI'])
        ad = ray.all_data()
        print((ad[('gas', 'O_p5_number_density')] * ad[('gas', 'dl')]).sum())

    Depending on the ion, you usually need to see at least :math:`10^{12} cm^{-2}`
    to have any appreciable absorption.  Try sending a sightline through a
    denser region in your simulation that might have more of that ion.

I don't have a metallicity field in my dataset--What can I do?
--------------------------------------------------------------

In order to estimate the density of ions throughout your dataset, Trident
needs a metallicity field.  But some datasets may not have one generated
by default.  I highly recommend re-running the dataset with metals present,
as this will lead to the best estimate of ions from Trident, but if you just
want to create a "dummy" metallicity field, include the following code at the
top of your script to automatically add a uniform metallicity field to any
datasets loaded lacking one (in this case it's 0.3 solar metallicity).  For more
information on creating derived fields like this one, see the `yt documentation
on derived fields
<https://yt-project.org/docs/dev/developing/creating_derived_fields.html>`_
::

        import yt
        import numpy as np

        def _metallicity(field, data):
            factor = 0.3 # 0.3 solar metallicity
            return (
                data.ds.arr(np.ones_like(data["gas", "density"]), 'Zsun') * factor
            )

        yt.add_field(
            ("gas", "metallicity"),
            function=_metallicity,
            sampling_type="local",
            units="Zsun",
        )

What functions are available and what is their syntax?
------------------------------------------------------

Go see the full documentation for all of our available classes and functions in the
:ref:`API Documentation <api-reference>`.

What version of Trident am I running?
-------------------------------------

To learn what version of Trident you're running, type::

    $ python
    >>> import trident
    >>> print(trident.__version__)

If you have a version ending in dev, it means you're on the development branch
and you should also figure out which particular changeset you're running.  You
can do this by::

    $ cd <path/to/trident>
    $ git log --pretty=format:'%h' -n 1

To figure out what version of yt you're running, type::

    $ yt version

If you're writing to the mailing list with a problem, be sure to include all
of the above with your bug report or question.

.. _where-installed:

Where is Trident installed?  Where are its data files?
------------------------------------------------------

One can easily identify where Trident is installed::

    $ python
    >>> import trident
    >>> print(trident.path)

The data files are located in that path with an appended ``/data``.

.. _mailing-list:

How do I join the mailing list?
-------------------------------

You can join our mailing list for announcements, bugs reports, and changes
at:

https://groups.google.com/forum/#!forum/trident-project-users

How do I learn more about the algorithms used in Trident?
---------------------------------------------------------

We have a full description of all the methods used in Trident including
citations to previous related works in our `Trident method paper
<http://adsabs.harvard.edu/abs/2017ApJ...847...59H>`_.

How do I cite Trident in my research?
-------------------------------------

Check out our :ref:`citation <citation>` page.

.. _slack-channel:

How do I get an invite to the Trident slack channel?
----------------------------------------------------

Click on this `link <https://join.slack.com/t/trident-project/shared_invite/zt-42h0uuwy-fBggZbeymnq2cB9ivtWloA>`_.
