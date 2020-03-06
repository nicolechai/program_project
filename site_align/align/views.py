
# Create your views here.

from django.http import HttpResponse, Http404, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic

from .models import Question

class IndexView(generic.ListView):
	model = Question
	template_name = './align/index.html'
	def index(request):
		question_list = Question.objects.all()
		output = ', '.join([q.question_text for q in question_list])
		return HttpResponse(output)


class DetailView(generic.DetailView):
    model = Question
    template_name = './align/mycgene.html'
