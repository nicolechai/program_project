from django.db import models

# Create your models here.

import datetime

from django.db import models
from django.utils import timezone

class Question(models.Model):
	question_text = models.CharField(max_length=300)
	def __str__(self):
		return self.question_text
