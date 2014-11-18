clear model am;

model.id = 'Bcl2_reversible'; % Bcl2 reversible model with Bak multimer and 
model.species = struct('id', {'B', 'T', 'M', 'D', 'aB', 'aT', 'aM', 'aB4'});
model.reaction(1) = struct('id', 'T+B->T+aB', 'reactant', struct('species', {'T' 'B'}, 'stoichiometry', {1 1}), 'product', struct('species', {'T' 'aB'}, 'stoichiometry', {1 1}), 'reversible', false);
model.reaction(2) = struct('id', 'aB+M<->aM', 'reactant', struct('species', {'aB' 'M'}, 'stoichiometry', {1 1}), 'product', struct('species', 'aM', 'stoichiometry', 1), 'reversible', true);
model.reaction(3) = struct('id', 'T+M<->aT', 'reactant', struct('species', {'T' 'M'}, 'stoichiometry', {1 1}), 'product', struct('species', 'aT', 'stoichiometry', 1), 'reversible', true);
model.reaction(4) = struct('id', '2aB<->D', 'reactant', struct('species', {'aB'}, 'stoichiometry', {2}), 'product', struct('species', 'D', 'stoichiometry', 1), 'reversible', true);
model.reaction(5) = struct('id', 'B+aB->2aB', 'reactant', struct('species', {'aB' 'B'}, 'stoichiometry', {1 1}), 'product', struct('species', 'aB', 'stoichiometry', 2), 'reversible', false);
model.reaction(5) = struct('id', '2D<->aB4', 'reactant', struct('species', {'D'}, 'stoichiometry', {2}), 'product', struct('species', {'aB4'}, 'stoichiometry', {1}), 'reversible', true);
model.reaction(6) = struct('id', 'aB->B', 'reactant', struct('species', {'aB'}, 'stoichiometry', {1}), 'product', struct('species', {'B'}, 'stoichiometry', {1}), 'reversible', false);
model.reaction(7) = struct('id', '0<->B', 'reactant', struct([]), 'product', struct('species', {'B'}, 'stoichiometry', {1}), 'reversible', 'true');
model.reaction(8) = struct('id', '0<->T', 'reactant', struct([]), 'product', struct('species', {'T'}, 'stoichiometry', {1}), 'reversible', 'true');
model.reaction(9) = struct('id', '0<->M', 'reactant', struct([]), 'product', struct('species', {'M'}, 'stoichiometry', {1}), 'reversible', 'true');
model.reaction(10) = struct('id', 'D->0', 'reactant', struct('species', {'D'}, 'stoichiometry', {1}), 'product', struct([]), 'reversible', 'false');
model.reaction(11) = struct('id', 'aB->0', 'reactant', struct('species', {'aB'}, 'stoichiometry', {1}), 'product', struct([]), 'reversible', 'false');
model.reaction(12) = struct('id', 'aM->0', 'reactant', struct('species', {'aT'}, 'stoichiometry', {1}), 'product', struct([]), 'reversible', 'false');
model.reaction(13) = struct('id', 'aT->0', 'reactant', struct('species', {'aM'}, 'stoichiometry', {1}), 'product', struct([]), 'reversible', 'false');
model.reaction(14) = struct('id', 'aB4->0', 'reactant', struct('species', {'aB4'}, 'stoichiometry', {1}), 'product', struct([]), 'reversible', 'false');
am = model_analysis(model);
