from flask import Blueprint, request, render_template, redirect, url_for
from flask_security import login_required
from lbrc_flask.database import db
from lbrc_flask.forms import SearchForm
from acadex.model import Academic
from scholarly import scholarly


blueprint = Blueprint("ui", __name__, template_folder="templates")

# Login required for all views
@blueprint.before_request
# @login_required
def before_request():
    pass


@blueprint.record
def record(state):
    if db is None:
        raise Exception(
            "This blueprint expects you to provide " "database access through database"
        )

@blueprint.route("/")
def index():
    search_form = SearchForm(formdata=request.args)

    q = Academic.query

    if search_form.search.data:
        q = q.filter(Academic.name.like("%{}%".format(search_form.search.data)))

    q = q.order_by(Academic.name.asc())

    academics = q.paginate(
            page=search_form.page.data,
            per_page=5,
            error_out=False,
        )

    return render_template("ui/index.html", academics=academics, search_form=search_form)


@blueprint.route("/add_search")
def add_search():
    search_form = SearchForm(formdata=request.args)

    academics = []

    if search_form.search.data:
        for a in scholarly.search_author(search_form.search.data):
            academics.append(a)

    return render_template("ui/add.html", academics=academics, search_form=search_form)


@blueprint.route("/add/<string:google_scholar_id>/")
def add(google_scholar_id):
    resp = scholarly.fill(scholarly.search_author_id(google_scholar_id), sections=['indices'])

    if resp:
        db.session.add(
            Academic(
                google_scholar_id=google_scholar_id,
                name=resp['name'],
                affiliation=resp['affiliation'],
                cited_by=resp['citedby'],
                h_index=resp['hindex'],
                i10_index=resp['i10index'],
        ))

        db.session.commit()

    return redirect(url_for('ui.index'))
