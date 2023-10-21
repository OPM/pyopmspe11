-- Copyright (C) 2023 NORCE

COORD
% for i in range(dic['noCells'][0] + 1):
 ${f"{dic['xcor'][i * (dic['noCells'][2] + 1)] : E}"} 0 ${f"{dic['zcor'][i*(dic['noCells'][2]+1)] : E}"}
 ${f"{dic['xcor'][i * (dic['noCells'][2] + 1) + dic['noCells'][2]] : E}"} 0 ${f"{dic['zcor'][((i + 1) * (dic['noCells'][2] + 1)) - 1] : E}"}
% endfor
% for i in range(dic['noCells'][0] + 1):
 ${f"{dic['xcor'][i * (dic['noCells'][2] + 1)] : E}"} ${dic['dims'][1]} ${f"{dic['zcor'][i*(dic['noCells'][2]+1)] : E}"}
 ${f"{dic['xcor'][i * (dic['noCells'][2] + 1) + dic['noCells'][2]] : E}"} ${dic['dims'][1]} ${f"{dic['zcor'][((i + 1) * (dic['noCells'][2] + 1)) - 1] : E}"}
% endfor
/

ZCORN
% for j in range(2):
% for i in range(dic['noCells'][0]):
 ${f"{dic['zcor'][i*(dic['noCells'][2]+1)] : E}"} ${f"{dic['zcor'][(i+1)*(dic['noCells'][2]+1)] : E}"}
% endfor
% endfor
% for k in range(dic['noCells'][2] - 1):
% for h in range(2):
% for j in range(2):
% for i in range(dic['noCells'][0]):
 ${f"{dic['zcor'][(i * (dic['noCells'][2] + 1)) + k + 1] : E}"} ${f"{dic['zcor'][((i + 1) * (dic['noCells'][2] + 1)) + k + 1] : E}"}
% endfor
% endfor
% endfor
% endfor
% for j in range(2):
% for i in range(dic['noCells'][0]):
 ${f"{dic['zcor'][((i + 1) * (dic['noCells'][2] + 1)) - 1] : E}"} ${f"{dic['zcor'][((i + 2) * (dic['noCells'][2] + 1)) - 1] : E}"}
% endfor
% endfor
/