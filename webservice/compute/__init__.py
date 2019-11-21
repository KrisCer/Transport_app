import random
import logging
import flask
from flask import Blueprint
from flask import request
import matplotlib.pyplot as plt
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import os
import time
directory = os.path.abspath(os.path.split(os.path.realpath(__file__))[0] + "/../")
template_folder = os.path.join(directory, 'templates/user_templates')
static_folder = os.path.join(directory, "static")
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
    energies,trans,trans2=plot_conductance(fsys,perfsys,emin=emin,emax=emax)
    plt.figure(figsize=(10,10))
    plt.xlabel("Energy (t)")
    plt.ylabel("Conductance (e\N{SUPERSCRIPT TWO}/h)")
    plt.xticks(np.arange(emin, emax+(emax-emin)/10, round((emax-emin)/10,2)),rotation=90)
    plt.grid(b=None, which='major', axis='both')
    plt.plot(energies, trans,label='Junction')
    plt.plot(energies, trans2,label='GNR')
    plt.legend()
    plt.savefig('webservice/user_static/img/cond.png',bbox_inches='tight')
    energies,dos,dos2=get_DOS(fsys,perfsys,emin=emin,emax=emax)
    plt.figure(figsize=(10,10))
    plt.plot(energies, dos,label='Junction')
    plt.plot(energies, dos2,label='GNR')
    plt.xlabel("Energy (t)")
    plt.ylabel("DOS")
    plt.xticks(np.arange(emin, emax+(emax-emin)/10, round((emax-emin)/10,2)),rotation=90)
    plt.grid(b=None, which='major', axis='both')
    plt.legend()
    plt.savefig('webservice/user_static/img/DOS.png',bbox_inches='tight')
    calculate='/user_static/img/Calculate_add.png'
    return flask.render_template('user_templates/app_page.html',band=calculate ,current=calculate,ldos=calculate,wfn=calculate,structure ='/user_static/img/plot.png'+'?'+str(time.time()), DOS='/user_static/img/DOS.png'+'?'+str(time.time()) ,cond ='/user_static/img/cond.png'+'?'+str(time.time()))

@blueprint.route("/compute/single_en_prop/",methods=["POST","GET"])
def display_band():
    momenta,energies=get_band(lead0)
    plt.figure(figsize=(10,10))
    plt.plot(momenta, energies)
    plt.xticks(np.arange(-1,1.1,1), ('X', 'Γ', 'X'))
    plt.xlabel("Wave vector")
    plt.ylabel("Energy (t)")
    plt.grid(b=None, which='major', axis='both')
    plt.savefig('webservice/user_static/img/band.png',bbox_inches='tight')
    en=float(request.form['en'])
    get_ldos(fsys,xaxis,yaxis,en)
    mode=int(request.form['mode'])
    lead_wfn(perfsys,en,mode)
    get_current(fsys,xaxis,yaxis,en)
    return flask.render_template('user_templates/app_page.html',current ='/user_static/img/current.png'+'?'+str(time.time()), wfn ='/user_static/img/wfn.png'+'?'+str(time.time()), structure ='/user_static/img/plot.png'+'?'+str(time.time()), DOS='/user_static/img/DOS.png'+'?'+str(time.time()) ,cond ='/user_static/img/cond.png'+'?'+str(time.time()), band='/user_static/img/band.png'+'?'+str(time.time()),ldos ='/user_static/img/ldos.png'+'?'+str(time.time()) )


