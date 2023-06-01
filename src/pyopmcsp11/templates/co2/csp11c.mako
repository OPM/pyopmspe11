-- Copyright (C) 2023 NORCE

----------------------------------------------------------------------------
RUNSPEC
----------------------------------------------------------------------------
DIMENS 
${dic['noCells'][0]} ${dic['noCells'][1]} ${dic['noCells'][2]} /

EQLDIMS
/

TABDIMS
${dic['noSands']} 1* ${dic['tabdims']} ${dic['tabdims']} /

OIL
GAS
CO2STORE
% if dic['model'] == 'complete':
DISGAS
DIFFUSE
THERMAL 
% endif

METRIC

START
1 'JAN' 2025 /

WELLDIMS
${len(dic['wellijk'])} ${1+max(dic['wellijkf'][0][1]-dic['wellijk'][0][1], dic['wellijkf'][1][1]-dic['wellijk'][1][1])} ${len(dic['wellijk'])} ${len(dic['wellijk'])} /

UNIFIN
UNIFOUT
----------------------------------------------------------------------------
GRID
----------------------------------------------------------------------------
INIT

INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/GRID.INC' /

INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/PERMX.INC' /

COPY 
PERMX PERMY /
PERMX PERMZ /
/

OPERATE
PERMZ 1* 1* 1* 1* 1* 1* 'MULTX' PERMZ 0.1 /
/

INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/PORO.INC' /

% if dic['model'] == 'complete':
INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/THCONR.INC' /
% endif
----------------------------------------------------------------------------
PROPS
----------------------------------------------------------------------------
INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/TABLES.INC' /

% if dic['model'] == 'complete':
DIFFC
18.01528 44.01 ${dic["diffusion"][1]} 1* ${dic["diffusion"][0]}  /

SPECROCK
% for i in range(dic['noSands']): 
${dic["temperature"][1]} ${dic["rockExtra"][0]}
${dic["temperature"][0]} ${dic["rockExtra"][0]} /
% endfor
% endif
----------------------------------------------------------------------------
REGIONS
----------------------------------------------------------------------------
INCLUDE
'${dic['exe']}/${dic['fol']}/preprocessing/SATNUM.INC' /
----------------------------------------------------------------------------
SOLUTION
---------------------------------------------------------------------------
EQUIL
 0 ${dic['pressure']/1.E5} ${dic['dims'][2] + dic['elevation']} 0 0 0 1 1 0 /

RPTRST
% if dic['model'] == 'immiscible': 
 'BASIC=2' FLOWS FLORES DEN/
% else:
 'BASIC=2' DEN/
%endif

% if dic['model'] == 'complete':
RSVD
0   0.0
${dic['dims'][2] + dic['elevation']} 0.0 /

RVVD
0   0.0
${dic['dims'][2] + dic['elevation']} 0.0 /

RTEMPVD
0   ${dic["temperature"][1]}
${dic['dims'][2] + dic['elevation']} ${dic["temperature"][0]} /
% endif
----------------------------------------------------------------------------
SCHEDULE
----------------------------------------------------------------------------
RPTRST
% if dic['model'] == 'immiscible': 
 'BASIC=2' FLOWS FLORES DEN/
% else:
 'BASIC=2' DEN/
%endif

WELSPECS
% for i in range(len(dic['wellijk'])):
'INJ${i}'	'G1'	${dic['wellijk'][i][0]}	${dic['wellijk'][i][1]}	1*	'GAS' ${dic['radius'][i]}/
% endfor
/
COMPDAT
% for i in range(len(dic['wellijk'])):
% for j in range(1+dic['wellijkf'][i][1]-dic['wellijk'][i][1]):
% if i==0:
'INJ${i}'	${dic['wellijk'][i][0]}	${dic['wellijk'][i][1]+j} ${dic['wellijk'][i][2]}	${dic['wellijk'][i][2]}	'OPEN'	1*	1*	${2.*dic['radius'][i]} /
% else:
'INJ${i}'	${dic['wellijk'][i][0]} ${dic['wellijk'][i][1]+j}	${dic['wellkh'][j]} ${dic['wellkh'][j]}	'OPEN'	1*	1*	${2.*dic['radius'][i]} /
%endif
% endfor
% endfor
/

% for j in range(len(dic['inj'])):
TUNING
1e-2 ${dic['inj'][j][2] / 86400.} 1e-10 2* 1e-12/
/
/
WCONINJE
% for i in range(len(dic['wellijk'])):
% if dic['inj'][j][3+3*i] > 0:
'INJ${i}' 'GAS' ${'OPEN' if dic['inj'][j][4+3*i] > 0 else 'SHUT'}
'RATE' ${f"{dic['inj'][j][4+3*i] * 86400 / 1.86843 : E}"}  1* 400/
% else:
'INJ${i}' 'OIL' ${'OPEN' if dic['inj'][j][4+3*i] > 0 else 'SHUT'}
'RATE' ${f"{dic['inj'][j][4+3*i] * 86400 / 998.108 : E}"}  1* 400/
%endif
% endfor
/
% if dic['model'] == 'complete':
WTEMP
% for i in range(len(dic['wellijk'])):
 'INJ${i}' ${dic['inj'][j][5+3*i]} /
% endfor
/
%endif
TSTEP
${round(dic['inj'][j][0]/dic['inj'][j][1])}*${dic['inj'][j][1] / 86400.}
/
% endfor