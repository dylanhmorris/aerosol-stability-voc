#####################################
# name: Makefile
# author: Dylan Morris <dylan@dylanhmorris.com>
#
# Makefile to generate analyses
# for project
####################################

#####################################
# Directory structure
####################################

default: all

SRC = src
OUTPUT = out
DATA = dat

RAW = $(DATA)/raw
CLEANED = $(DATA)/cleaned
PARAMS = $(SRC)/parameters
STAN = $(SRC)/stan

CHAINS = $(OUTPUT)/chains
FIGDIR = $(OUTPUT)/figures
TABDIR = $(OUTPUT)/tables

#####################################
# Expected shell settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

MKDIR := @mkdir -p
RM := rm -f
RMDIR := rmdir
CP := @cp
LATEXMK := latexmk -xelatex

HIDDEN := .DS_Store

R_OPTIONS = --no-save --no-restore --no-site-file --no-environ
R_COMMAND := Rscript $(R_OPTIONS)

# R creates a blank Rplots.pdf when run
# from the command line. This removes it.
FIG_CLEANUP = @$(RM) Rplots.pdf

#####################################
# Code/data locations
#####################################
## enumerate data
DATAFILE = data
PCR_FILE = pcr
LINEAGE_DATAFILE = pango_substitutions
SURFACE_DATAFILE = surface-data
DATA_EXT = .tsv

DATAFILES = $(DATAFILE) $(SURFACE_DATAFILE) $(LINEAGE_DATAFILE)

DATA_PATHS = $(addsuffix $(DATA_EXT), $(addprefix $(CLEANED)/, $(DATAFILES)))

HYPERS = $(PARAMS)/hypers.R

TITER_MODEL = $(STAN)/infer_titers.stan
AEROSOL_MODEL = $(STAN)/infer_timeseries_halflives.stan
SURFACE_MODEL = $(STAN)/infer_halflives.stan


## enumerate data cleaning scripts
AEROSOL_CLEANING_SCRIPT = $(SRC)/clean-aerosol-data.R
SURFACE_CLEANING_SCRIPT = $(SRC)/clean-surface-data.R

MODEL_NAMES = titers-aerosol halflives-aerosol halflives-aerosol-no-pcr titers-surface halflives-surface

## enumerate all figures
PRIOR_CHECK_FIGS = figure-prior-check-aerosol.pdf figure-prior-check-surface.pdf

POSTERIOR_CHECK_FIGS = figure-posterior-check-aerosol-vero.pdf figure-posterior-check-aerosol-tmprss2-reg.pdf figure-posterior-check-aerosol-tmprss2-rml.pdf figure-posterior-check-surface.pdf

MAIN_ANALYSIS_FIGS = figure-main.pdf figure-main-surface.pdf figure-aa-subst.pdf

SENSITIVITY_FIGS = figure-main-no-pcr.pdf  figure-main-linear-surface.pdf figure-tmprss2-reg.pdf figure-tmprss2-rml.pdf figure-no-pcr-tmprss2-reg.pdf figure-no-pcr-tmprss2-rml.pdf

CELL_LINE_FIGS = figure-halflives.pdf figure-halflives-no-pcr.pdf figure-ratios.pdf figure-ratios-no-pcr.pdf

RAW_TIMESERIES_FIGS = figure-raw-timeseries-vero.pdf figure-raw-timeseries-tmprss2-reg.pdf figure-raw-timeseries-tmprss2-rml.pdf

FIGURE_NAMES = $(PRIOR_CHECK_FIGS) $(POSTERIOR_CHECK_FIGS) $(MAIN_ANALYSIS_FIGS) $(SENSITIVITY_FIGS) $(CELL_LINE_FIGS) $(RAW_TIMESERIES_FIGS)


TABLE_NAMES = table-halflives-aerosol.tsv table-halflives-surface.tsv table-halflives-aerosol-no-pcr.tsv table-ratios-aerosol.tsv table-ratios-surface.tsv table-ratios-aerosol-no-pcr.tsv

CHECK_PATHS = $(addprefix $(CHAINS)/, $(addsuffix -prior-check.Rds, $(MODEL_NAMES)))
CHAIN_PATHS = $(addprefix $(CHAINS)/, $(addsuffix -chains.Rds, $(MODEL_NAMES)))
PRIOR_CHECK_FIGURE_PATHS = $(addprefix $(FIGDIR)/, $(PRIOR_CHECK_FIGS))
FIGURE_PATHS = $(addprefix $(FIGDIR)/, $(FIGURE_NAMES))
TABLE_PATHS = $(addprefix $(TABDIR)/, $(TABLE_NAMES))

#####################################
# Rules
#
# definition of dependency
# tree and specification of
# rules for doing stuff
#####################################

##########################
# rules for data cleaning
##########################

$(CLEANED)/$(LINEAGE_DATAFILE)$(DATA_EXT): $(RAW)/$(LINEAGE_DATAFILE)$(DATA_EXT)
	$(MKDIR) $(CLEANED)
	$(CP) $< $@

$(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT): $(SURFACE_CLEANING_SCRIPT) $(RAW)/$(SURFACE_DATAFILE)$(DATA_EXT)
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@


$(CLEANED)/%$(DATA_EXT): $(AEROSOL_CLEANING_SCRIPT) $(RAW)/$(DATAFILE).csv $(RAW)/$(PCR_FILE).csv
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@

##########################
# rules for fitting
##########################

$(CHAINS)/halflives-aerosol-no-pcr-%.Rds: $(SRC)/infer-aerosol-halflives.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(HYPERS) $(AEROSOL_MODEL)
	$(MKDIR) $(CHAINS)
	$(R_COMMAND) $^ $@

$(CHAINS)/titers-surface-%.Rds: $(SRC)/infer-titers.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(HYPERS)  $(TITER_MODEL)
	$(MKDIR) $(CHAINS)
	$(R_COMMAND) $^ $@

$(CHAINS)/halflives-surface-%.Rds: $(SRC)/infer-surface-halflives.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(HYPERS) $(SURFACE_MODEL)
	$(MKDIR) $(CHAINS)
	$(R_COMMAND) $^ $@

$(CHAINS)/halflives-aerosol-%.Rds: $(SRC)/infer-aerosol-halflives.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(HYPERS) $(AEROSOL_MODEL)
	$(MKDIR) $(CHAINS)
	$(R_COMMAND) $^ $@

$(CHAINS)/titers-aerosol-%.Rds: $(SRC)/infer-titers.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(HYPERS) $(TITER_MODEL)
	$(MKDIR) $(CHAINS)
	$(R_COMMAND) $^ $@


##########################
# rules for figures
##########################
$(FIGDIR)/figure-aa-subst.pdf: $(SRC)/figure-aa-subst.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-chains.Rds $(CLEANED)/$(LINEAGE_DATAFILE).tsv
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-posterior-check-aerosol-%.pdf: $(SRC)/figure-predictive-check.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/titers-aerosol-chains.Rds $(CHAINS)/halflives-aerosol-chains.Rds 
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-prior-check-aerosol.pdf: $(SRC)/figure-predictive-check.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/titers-aerosol-chains.Rds $(CHAINS)/halflives-aerosol-prior-check.Rds 
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-posterior-check-surface.pdf: $(SRC)/figure-predictive-check.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/titers-surface-chains.Rds $(CHAINS)/halflives-surface-chains.Rds 
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-prior-check-surface.pdf: $(SRC)/figure-predictive-check.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/titers-surface-chains.Rds $(CHAINS)/halflives-surface-prior-check.Rds 
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-no-pcr-tmprss2-%.pdf: $(SRC)/figure-main.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-no-pcr-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-tmprss2-%.pdf: $(SRC)/figure-main.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/%-no-pcr.pdf: $(SRC)/%.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-no-pcr-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/figure-main-linear-surface.pdf: $(SRC)/figure-main.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-surface-chains.Rds $(CHAINS)/titers-surface-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGDIR)/%-surface.pdf: $(SRC)/%.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-surface-chains.Rds $(CHAINS)/titers-surface-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(FIGDIR)/figure-raw-timeseries-%.pdf: $(SRC)/figure-raw-timeseries.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(FIGDIR)/%.pdf: $(SRC)/%.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(FIGDIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


##########################
# rules for tables
##########################
$(TABDIR)/table-halflives-aerosol.tsv: $(SRC)/figure-halflives.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@

$(TABDIR)/table-halflives-surface.tsv: $(SRC)/figure-halflives.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-surface-chains.Rds $(CHAINS)/titers-surface-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@

$(TABDIR)/table-halflives-aerosol-no-pcr.tsv: $(SRC)/figure-halflives.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-no-pcr-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@

$(TABDIR)/table-ratios-aerosol.tsv: $(SRC)/figure-ratios.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@


$(TABDIR)/table-ratios-surface.tsv: $(SRC)/figure-ratios.R $(CLEANED)/$(SURFACE_DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-surface-chains.Rds $(CHAINS)/titers-surface-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@


$(TABDIR)/table-ratios-aerosol-no-pcr.tsv: $(SRC)/figure-ratios.R $(CLEANED)/$(DATAFILE)$(DATA_EXT) $(CHAINS)/halflives-aerosol-no-pcr-chains.Rds $(CHAINS)/titers-aerosol-chains.Rds
	$(MKDIR) $(TABDIR)
	$(R_COMMAND) $^ $@


# various convenience commands
.PHONY: deltemp clean install uninstall renv destroyrenv rebuild freeze data figures chains priorchecks all

## remove emacs tempfiles, etc.
deltemp:
	$(RM) $(SRC)/*~*
	$(RM) $(SRC)/*#*
	$(RM) $(PARAMS)/*~*
	$(RM) $(PARAMS)/*#*

clean: deltemp
	$(RM) $(CLEANED)/*.tsv $(CLEANED)/$(HIDDEN)
	$(MKDIR) $(CLEANED)
	$(RMDIR) $(CLEANED)
	$(RM) $(CHAINS)/*.Rds $(CHAINS)/$(HIDDEN)
	$(RM) $(FIGDIR)/*.pdf $(FIGDIR)/*.png $(FIGDIR)/$(HIDDEN)
	$(RM) $(TABDIR)/*.tsv $(TABDIR)/$(HIDDEN)
	$(MKDIR) $(CHAINS) $(FIGDIR) $(TABDIR)
	$(RMDIR) $(CHAINS)
	$(RMDIR) $(FIGDIR)
	$(RMDIR) $(TABDIR)
	$(RM) $(OUTPUT)/$(HIDDEN)
	$(MKDIR) $(OUTPUT)
	$(RMDIR) $(OUTPUT)
	$(RM) $(STAN)/*.rds $(STAN)/$(HIDDEN)


install: $(SRC)/install.R
	@echo "\nAttempting to install project package and dependecies...\n\n"
	$(R_COMMAND) $^


uninstall: $(SRC)/uninstall.R
	@echo "\nAttempting to uninstall project package...\n\n"
	$(R_COMMAND) $^

renv: destroyrenv
	$(MKDIR) renv/cellar
	$(R_COMMAND) -e "devtools::build('varstab', path='renv/cellar')"
	$(R_COMMAND) -e "renv::init()"
	$(R_COMMAND) -e "remotes::install_local('varstab', force=TRUE)"
	$(R_COMMAND) -e "renv::snapshot()"
	@echo '\ncat("Virtual environment activated successfully for project *variant aerosol stability*.\n")\n' >> .Rprofile

freeze: $(SRC)/collect-dependencies.R
	$(R_COMMAND) $<

destroyrenv:
	@echo "\nAttempting to destroy virtual environment...\n\n"
	$(RM) .Rprofile
	$(RM) renv.lock
	$(RM) -r renv/*
	$(RM) renv/.gitignore
	$(RM) renv/.DS_Store
	$(MKDIR) renv
	$(RMDIR) renv

rebuild: clean uninstall install all

## group target aliases
data: $(DATA_PATHS)
figures: $(FIGURE_PATHS)
tables: $(TABLE_PATHS)
chains: $(CHAIN_PATHS) $(CHECK_PATHS) 
priorchecks: $(PRIOR_CHECK_FIGURE_PATHS)
all: $(DATA_PATHS) $(CHECK_PATHS) $(CHAIN_PATHS) $(FIGURE_PATHS) $(TABLE_PATHS)
