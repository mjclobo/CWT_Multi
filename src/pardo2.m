function out1 = pardo2(outdim1,fcn)

out1 = cell(outdim1,1);

it = 1:outdim1;

for i = it
    out1{i} = fcn(i);
end

end
