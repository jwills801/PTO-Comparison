function [P_Map Q_Map W_Map Disp Cf Ch Cv Cs Cst] = PumpInfo(PumpName)

eval(strjoin(['load '  string(PumpName)]))

end