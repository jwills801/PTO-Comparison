clear
n = 3;

% Make Symbolic variables
syms t 
symstr = 'syms ';
Mstr = 'syms';
Cstr = 'syms';
Kstr = 'syms';
Fexcstr = 'syms';
Fptostr = 'syms';
for i = 1:n
        symstr = [symstr,'x',num2str(i),'(t) '];
        Fexcstr = [Fexcstr,' Fexc',num2str(i)];
        Fptostr = [Fptostr,' Fpto',num2str(i)];
        for j = 1:n
            Mstr = [Mstr , ' M',num2str(i),num2str(j)];
            Cstr = [Cstr , ' C',num2str(i),num2str(j)];
            Kstr = [Kstr , ' K',num2str(i),num2str(j)];
        end
end
eval(symstr), eval(Mstr), eval(Cstr), eval(Kstr), eval(Fexcstr), eval(Fptostr)

% Put into matrices
for i = 1:n
    x(i,1) = eval(['x',num2str(i)]);
    Fexc(i,1) = eval(['Fexc',num2str(i)]);
    Fpto(i,1) = eval(['Fpto',num2str(i)]);
    for j = 1:n
        M(i,j) = eval(['M',num2str(i),num2str(j)]);
        C(i,j) = eval(['C',num2str(i),num2str(j)]);
        K(i,j) = eval(['K',num2str(i),num2str(j)]);
    end
end

% Equations of motion:
EOM = M*diff(x,t,t) + C*diff(x,t) + K*x - Fexc - Fpto == zeros(size(x))


