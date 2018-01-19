#        Lower-level Makefile for RIP

include ../configure.rip

OBJS= \
	addorfill.o adjellip.o angle.o barb.o bes.o bogs.o \
	brnshr.o bshear.o bvfricalc.o capecalc3d.o \
	ceilingcalc.o colram.o colrln.o condheat.o contrive2.o \
	cpchhl.o cpchll.o \
	cpcolr.o cpmpxy.o createdtg.o cttcalc.o dbzcalc.o refmos2calc.o refmos3calc.o ddpcalc.o \
	ddx.o ddy.o derivc.o \
	drawcl.o driver.o egridfill.o eqthecalc.o \
	extingcalc.o fields.o fillarray.o finterp.o \
	finterp2d.o fregcalc.o frtitle.o fx.o fy.o getbasicvars.o \
	getbasicvars_newdom.o getchar.o getconvals.o getdash.o getdims.o \
	getheadinfo.o getheadinfo_newdom.o getinum.o \
	getminfo.o getpt.o getrnum.o gettimes.o getvar.o getvarinfo.o \
	hboxdraw.o hbulldraw.o \
	hchadraw.o hcondraw.o hiresmap.o hlinedraw.o \
	hmapdraw.o hodograph.o hsidsdraw.o hstrdraw.o \
	hticdraw.o htrajdraw.o hvecdraw.o igetcoind.o lennonblank.o \
	liftparcel.o locinterp.o maptform.o maptick.o \
	mconvert.o mktrjpts.o mlen.o monotonic.o netasc.o \
	numeric.o omgcalc.o pfcalc.o pltitle.o prcver.o precprob.o \
	premaptform.o profil.o profvelcalc.o pvocalc.o qgomg.o rdcolt.o \
	readspec.o relhl.o rhucalc.o ricalc.o \
	rip.o rotipslab.o rotpslab.o sateqthecalc.o saweli.o scondraw.o seaprs.o \
        setripheader.o smooth.o smoothcp.o smoothcz.o smoothrain.o srflow.o \
	srflow2.o sstats.o sticdraw.o sticdraw_polar.o \
        svecdraw.o tdpcalc.o trop.o refalt.o lifted_index.o \
	tfpcalc.o thecalc.o tonpsadiabat.o \
	tserprep.o turb.o typhoon.o unblank.o vbardraw.o vc2dcalc.o vcondraw.o \
	velbrb.o velvctmts.o vgp.o vinterp.o \
	virtual.o virtualhyd.o viscalc.o vticdraw.o vtrajdraw.o \
	vvecdraw.o vwinddraw.o wdircalc.o wetbulbcalc.o \
	writefile.o wspcalc.o xmapcalc.o xtodot.o binio.o v5d.o rotate_earth2comp.o

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.

all: 	rip ripcomp ripdp_mm5 ripdp_wrfarw ripdp_mpas ripdp_wrfnmm ripcut ripinterp ripshow showtraj tabdiag upscale

rip: 	$(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(NETCDFLIB) $(NETCDFINC) $(LDFLAGS) $(NCARGLIBS) $(LOCAL_LIBS)

ripcomp:	ripcomp.f
	$(FC) $(FFLAGS) -o $@ $@.f 
ripdp_mm5:	ripdp_mm5.f	fillarray.o	mconvert.o	writefile_rdp.o	\
		virtual.o
	$(FC) $(FFLAGS) -o $@ $@.f fillarray.o mconvert.o writefile_rdp.o \
		virtual.o 
ripdp_wrfarw:	ripdp_wrfarw.f	commptf	premaptform.o	maptform.o	fillarray.o \
		mconvert.o	writefile_rdp.o	virtual.o	xtodot.o	commptf
	$(FC) $(FFLAGS) -o $@  $@.f $(NETCDFLIB)  $(NETCDFINC)  \
		premaptform.o maptform.o fillarray.o mconvert.o writefile_rdp.o \
		virtual.o xtodot.o 
ripdp_mpas:	ripdp_mpas.f	commptf premaptform.o   maptform.o      fillarray.o \
		mconvert.o	writefile_rdp.o virtual.o       xtodot.o        commptf
	$(FC) $(FFLAGS) -o $@  $@.f $(NETCDFLIB) $(NETCDFINC) \
		premaptform.o maptform.o fillarray.o mconvert.o writefile_rdp.o \
		virtual.o xtodot.o rotate_earth2comp.o
ripdp_wrfnmm:	ripdp_wrfnmm.f	commptf	maptform.o	premaptform.o	fillarray.o \
		mconvert.o	writefile_rdp.o	virtual.o	xmapcalc.o	egridfill.o	xtodot.o	commptf	CMASSI.comm
	$(FC) $(FFLAGS) -o $@  $@.f $(NETCDFLIB)  $(NETCDFINC) \
		maptform.o	premaptform.o fillarray.o mconvert.o writefile_rdp.o \
		virtual.o xmapcalc.o egridfill.o xtodot.o 
ripcut:	ripcut.f
	$(FC) $(FFLAGS) -o $@ $@.f 
ripinterp:	ripinterp.f
	$(FC) $(FFLAGS) -o $@ $@.f 
ripshow:	ripshow.f
	$(FC) $(FFLAGS) -o $@ $@.f 
showtraj:	showtraj.f
	$(FC) $(FFLAGS) -o $@ $@.f 
tabdiag:	tabdiag.f	lennonblank.f
	$(FC) $(FFLAGS) -o $@ $@.f lennonblank.f 
upscale:	upscale.f
	$(FC) $(FFLAGS) -o $@ $@.f 

addorfill.o:    comconst
adjellip.o:     comconst
binio.o:        binio.h
bogs.o:         comconst
brnshr.o:       comconst
bvfricalc.o:    comconst
capecalc3d.o:   comconst
ceilingcalc.o:  comconst
condheat.o:     comconst
contrive2.o:    comconst
cpmpxy.o:       comconst    comvctran
cttcalc.o:      comconst
dbzcalc.o:      comconst
refmos2calc.o:      comconst
	$(FC) $(FFLAGS) -c $(NETCDFLIB)  $(NETCDFINC) $*.f
refmos3calc.o:      comconst
	$(FC) $(FFLAGS) -c $(NETCDFLIB)  $(NETCDFINC) $*.f
ddpcalc.o:      comconst
ddx.o:          comconst
ddy.o:          comconst
derivc.o:       comconst
driver.o:       comconst    pointers commptf  v5df.h  getpt.o comvctran	driver.f
eqthecalc.o:    comconst
extingcalc.o:   comconst
fields.o:       comconst    pointers getpt.o fields.f
fregcalc.o:     comconst
fx.o:           comconst    comvctran
fy.o:           comconst    comvctran
getchar.o:      comconst
getdash.o:      comconst
getinum.o:      comconst
getpt.o:        comconst
getrnum.o:      comconst
getvar.o:       comconst
hboxdraw.o:     comconst
hbulldraw.o:    comconst
hchadraw.o:     comconst
hcondraw.o:     comconst	hcondraw.f
hiresmap.o:     comconst
hlinedraw.o:    comconst
hmapdraw.o:     comconst
hsidsdraw.o:    comconst
hstrdraw.o:     comconst
hticdraw.o:     comconst
htrajdraw.o:    comconst
hvecdraw.o:     comconst
igetcoind.o:    comconst
lennonblank.o:  comconst
lifted_index.o:   comconst
liftparcel.o:   comconst
locinterp.o:    comconst
maptform.o:	commptf
netasc.o:       comconst
omgcalc.o:      comconst
pfcalc.o:       comconst
pltitle.o:      comconst
prcver.o:       comconst
precprob.o:     comconst
premaptform.o:	commptf
profil.o:       comconst
profvelcalc.o:       comconst
pvocalc.o:      comconst
qgomg.o:        comconst	qgomg.f
rdcolt.o:       comconst
refalt.o: 
readspec.o:     comconst	readspec.f
rhucalc.o:      comconst
ricalc.o:      comconst
rip.o:          comconst
sateqthecalc.o:             comconst
saweli.o:       comconst
scondraw.o:     comconst
seaprs.o:
setripheader.o:	comconst
smooth.o:       comconst
smoothcp.o:     comconst
smoothcz.o:     comconst
sstats.o:       comconst	sstats.f
sticdraw.o:     comconst
sticdraw_polar.o:           comconst
svecdraw.o:     comconst
tdpcalc.o:      comconst
tfpcalc.o:      comconst
thecalc.o:      comconst
tonpsadiabat.o: comconst
trop.o:
tserprep.o:     comconst
turb.o:         comconst
typhoon.o:      comconst
unblank.o:      comconst
v5d.o:          binio.h     v5d.h       vis5d.h
vbardraw.o:	comconst
vc2dcalc.o:     comconst    comvctran
vcondraw.o:     comconst	vcondraw.f
velbrb.o:       comconst
vinterp.o:      comconst
virtual.o:
virtualhyd.o:
viscalc.o:      comconst
vticdraw.o:     comconst
vtrajdraw.o:    comconst
vvecdraw.o:     comconst    comvctran
vwinddraw.o:    comconst    comvctran
wdircalc.o:     comconst
wetbulbcalc.o:  comconst
writefile.o:    comconst
wspcalc.o:      comconst
xmapcalc.o:	commptf
rotate_earth2comp.o:	commptf

clean: 
	${RM} *.o
