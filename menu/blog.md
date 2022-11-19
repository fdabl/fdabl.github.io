---
layout: page
title: 
---

{% for post in site.posts %}
<div class="posts">
  <h1>
    <a href="{{ site.github.url | replace: 'http://', 'https://' }}{{ post.url }}">{{ post.title }}</a>
  </h1>
  {% if post.image %}
  <div class="thumbnail-container">
    <a href="{{ site.github.url | replace: 'http://', 'https://' }}{{ post.url }}"><img src="{{ site.github.url | replace: 'http://', 'https://' }}/assets/img/{{ post.image }}"></a>
  </div>
  {% endif %}
  <p>
    {{ post.content | strip_html | truncate: 350 }} <a href="{{ site.github.url | replace: 'http://', 'https://' }}{{ post.url }}">Read more</a>
    <span class="post-date"><i class="fa fa-calendar" aria-hidden="true"></i> {{ post.date | date_to_string }} - <i class="fa fa-clock-o" aria-hidden="true"></i> {% include read-time.html %}</span>
  </p>
</div>
{% endfor %}


<!-- Pagination links -->
<div class="pagination">
  {% if paginator.next_page %}
    <a class="pagination-button pagination-active next" href="{{ site.github.url | replace: 'http://', 'https://' }}}{{ paginator.next_page_path }}">{{ site.data.settings.pagination.previous_page }}</a>
  {% else %}
    <span class="pagination-button">{{ site.data.settings.pagination.previous_page }}</span>
  {% endif %}

  {% if paginator.previous_page %}
    <a class="pagination-button pagination-active" href="{{ site.baseurl| replace: 'http://', 'https://' }}}{{ paginator.previous_page_path }}">{{ site.data.settings.pagination.next_page }}</a>
  {% else %}
    <span class="pagination-button">{{ site.data.settings.pagination.next_page }}</span>
  {% endif %}
</div>
