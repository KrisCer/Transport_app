import random
import logging
import flask
from flask import Blueprint
from flask import request
import matplotlib.pyplot as plt
import io
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import os
import time
import kwant
import numpy as np
from tapp import make_system
from tapp import plot_numbers
from remove_atoms import del_atoms
from cond_DOS import plot_conductance
from cond_DOS import perfect_system
from cond_DOS import get_DOS
from ldos_current import get_band
from ldos_current import get_ldos
from ldos_current import lead_wfn
from ldos_current import get_current

blueprint = Blueprint('compute', __name__)

logger = logging.getLogger('tools-app')

@blueprint.route('/compute/show_system/', methods=["POST","GET"])
def return_fig():
    angle=int(request.form['angle'])
    width=int(request.form['width'])
    chiral=0
    if width%2==0:
        chiral=int(request.form['chiral'])
    shift=0
    if angle==120:
        distance=int(request.form['distance'])
    if angle==60:
        distance=1
    if angle==180:
        distance=int(request.form['distance'])
        shift=int(request.form['shift'])
    global sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W
    sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W = make_system(A=angle,W=width,D=distance,S=shift,F=chiral)
    global xaxis, yaxis
    xaxis, yaxis=plot_numbers(sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W)
    compute='/user_static/img/Compute_cond.png'
    return flask.render_template('user_templates/app_page.html',ldos=compute,current=compute,wfn=compute,band=compute,cond=compute,DOS=compute , structure ='/user_static/img/plot.png'+'?'+str(time.time()))

@blueprint.route('/compute/', methods=['POST','GET'])
def go_to_app_page():
    junction='/user_static/img/Design_junction.png'
    return flask.render_template('user_templates/app_page.html',band=junction, cond=junction ,DOS=junction, structure=junction,current=junction,ldos=junction,wfn=junction)

@blueprint.route("/compute/remove_atoms/",methods=["POST","GET"])
def return_fig_remove_atoms():
    atoms_to_remove=request.form['atoms_remove']
    atoms_to_remove=list(map(int, atoms_to_remove.split(',')))
    global sys
    sys=del_atoms(sys,a,b,atoms_to_remove)
    global xaxis, yaxis
    xaxis, yaxis=plot_numbers(sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W)
    compute='/user_static/img/Compute_cond.png'
    return flask.render_template('user_templates/app_page.html',ldos=compute,current=compute,wfn=compute,band=compute,cond=compute,DOS=compute,structure ='/user_static/img/plot.png'+'?'+str(time.time()))

@blueprint.route("/compute/conductivity_DOS/",methods=["POST","GET"])
def calculate_conductivity():
    emin=float(request.form['emin'])
    emax=float(request.form['emax'])
    global fsys
    fsys=sys.finalized()
    global perfsys
    perfsys=perfect_system(W)
    plot_conductance(fsys,perfsys,emin=emin,emax=emax)
    get_DOS(fsys,perfsys,emin=emin,emax=emax)
    calculate='/user_static/img/Calculate_add.png'
    return flask.render_template('user_templates/app_page.html',band=calculate ,current=calculate,ldos=calculate,wfn=calculate,structure ='/user_static/img/plot.png'+'?'+str(time.time()), DOS='/user_static/img/DOS.png'+'?'+str(time.time()) ,cond ='/user_static/img/cond.png'+'?'+str(time.time()))

@blueprint.route("/compute/single_en_prop/",methods=["POST","GET"])
def display_band():
    get_band(lead0) #Saves figure of band structure
    en=float(request.form['en'])
    get_ldos(fsys,xaxis,yaxis,en) #Saves figure of LDOS
    mode=int(request.form['mode'])
    lead_wfn(perfsys,en,mode) #Saves figure of wave function
    get_current(fsys,xaxis,yaxis,en) #Saves figure of current
    return flask.render_template('user_templates/app_page.html',current ='/user_static/img/current.png'+'?'+str(time.time()), wfn ='/user_static/img/wfn.png'+'?'+str(time.time()), structure ='/user_static/img/plot.png'+'?'+str(time.time()), DOS='/user_static/img/DOS.png'+'?'+str(time.time()) ,cond ='/user_static/img/cond.png'+'?'+str(time.time()), band='/user_static/img/band.png'+'?'+str(time.time()),ldos ='/user_static/img/ldos.png'+'?'+str(time.time()) )


