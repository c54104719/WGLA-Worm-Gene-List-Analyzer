from django.urls import path
from . import views

urlpatterns = [
    path('', views.home_view, name='index'),
    path('tool/', views.tool_view, name='tool'),
    path('tool2/', views.tool2_view, name='tool2'),
    path('download/', views.download_analysis_results, name='download_results'),
    path('help/', views.help_view, name='help'),
    path('developer/', views.developer_view, name='developer'),
]