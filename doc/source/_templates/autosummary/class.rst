{% extends "!autosummary/class.rst" %}

{% block methods %}
{% if methods %}

Methods
^^^^^^^

.. autosummary::
   :toctree:
   :nosignatures:
   {% for item in methods %}
       ~{{ name }}.{{ item }}
   {%- endfor %}

{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}

Attributes
^^^^^^^^^^

.. autosummary::
   :toctree:
   :nosignatures:
   {% for item in attributes %}
       ~{{ name }}.{{ item }}
   {%- endfor %}

{% endif %}
{% endblock %}
