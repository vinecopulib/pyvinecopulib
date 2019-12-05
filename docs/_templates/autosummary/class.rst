{{ objname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

   {% if objname != 'BicopFamily' %}
   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :nosignatures:
   {% for item in attributes %}
      ~{{ objname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :nosignatures:
      :toctree: .
   {% for item in methods %}
      ~{{ objname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
   {% endif %}
