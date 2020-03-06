from django.urls import path

from . import views

app_name = 'align'
urlpatterns = [
	# ex: /align/
	path('', views.IndexView.as_view(), name='index'),
	# ex: /align/5/
	path('<int:pk>/results', views.DetailView.as_view(), name='detail'),
]
