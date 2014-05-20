function model = compression(parameters)
	% import COMSOL features
	import com.comsol.model.*
	import com.comsol.model.util.*

	% create the model
	model = ModelUtil.create('Model');
	model.modelPath('compression');
	model.modelNode.create('mod1');
	model.geom.create('geom1', 2);
	model.mesh.create('mesh1', 'geom1');
	model.physics.create('solid', 'SolidMechanics', 'geom1');
	model.study.create('std1');
	model.study('std1').feature.create('stat', 'Stationary');
	model.study('std1').feature('stat').activate('solid', true);

	% set parameters in the model
	model.param.set('modelWidth', sprintf('%f[m]', parameters.domainWidth));
	model.param.set('modelDepth', sprintf('%f[m]', parameters.domainDepth));
	model.param.set('appliedStrain', sprintf('%f', parameters.appliedStrain));
	model.param.set('compression', 'modelDepth*appliedStrain');
	model.param.set('basalStiffness', sprintf('%f[Pa]', parameters.basalStiffness));
	model.param.set('stiffnessRatio', sprintf('%f', parameters.lesionStiffnessRatio));
	model.param.set('lesionStiffness', 'basalStiffness*stiffnessRatio');
	model.param.set('lesionExtraStiffness', 'lesionStiffness-basalStiffness');
	model.param.set('density', sprintf('%f[kg/m^3]', parameters.density));
	model.param.set('poissonsRatio', sprintf('%f', parameters.poissonsRatio));

	% deal with the human model case
	if strcmpi(cell2mat(parameters.caseCategory), 'human') == 1
		model.param.set('fatStiffness', '80[kPa]');
		model.param.set('boneStiffness', '18.6[GPa]');
		model.param.set('fatExtraStiffness', 'fatStiffness - basalStiffness');
		model.param.set('boneExtraStiffness', 'boneStiffness - basalStiffness');

		% import our stiffness image
		model.func.create('im1', 'Image');
		model.func('im1').set('filename', sprintf('compression/stiffnessMap_%03d.png', caseIndex));
		model.func('im1').importData;
		model.func('im1').set('xmin', sprintf('%f', (-parameters.domainWidth / 2)));
		model.func('im1').set('xmax', sprintf('%f', (parameters.domainWidth / 2)));
		model.func('im1').set('ymax', sprintf('%f', parameters.domainDepth));
		model.func('im1').set('inplace', 'off');
		model.func('im1').set('scaling', 'manual');
		model.func('im1').set('manualexpr', 'r');
		model.func('im1').set('funcname', 'lesionMap');
		model.func.duplicate('im2', 'im1');
		model.func('im2').set('funcname', 'boneMap');
		model.func('im2').set('manualexpr', 'b');
		model.func.duplicate('im3', 'im2');
		model.func('im3').set('funcname', 'fatMap');
		model.func('im3').set('manualexpr', 'g');
	else
		% import our stiffness image
		model.func.create('im1', 'Image');
		model.func('im1').set('funcname', 'stiffnessMap');
		model.func('im1').set('filename', sprintf('compression/stiffnessMap_%03d.png', caseIndex));
		model.func('im1').importData;
		model.func('im1').set('xmin', sprintf('%f', (-parameters.domainWidth / 2)));
		model.func('im1').set('xmax', sprintf('%f', (parameters.domainWidth / 2)));
		model.func('im1').set('ymax', sprintf('%f', parameters.domainDepth));
	end

	% define the model geometry
	model.geom('geom1').feature.create('r1', 'Rectangle');
	model.geom('geom1').feature('r1').setIndex('size', 'domainWidth', 0);
	model.geom('geom1').feature('r1').setIndex('size', 'domainDepth', 1);
	model.geom('geom1').feature('r1').setIndex('size', 'modelWidth', 0);
	model.geom('geom1').feature('r1').setIndex('size', 'modelDepth', 1);
	model.geom('geom1').feature('r1').setIndex('pos', '-modelWidth/2', 0);
	model.geom('geom1').run;

	% set the material properties
	model.physics('solid').feature('lemm1').set('NearlyIncompressible', 1, '1');
	model.physics('solid').feature('lemm1').set('E_mat', 1, 'userdef');

	% deal with the human model case
	if strcmpi(cell2mat(parameters.caseCategory), 'human') == 1
		model.physics('solid').feature('lemm1').set('E', 1, 'basalStiffness + (fatMap(x, modelDepth - y) * fatExtraStiffness) + (lesionMap(x, modelDepth - y) * lesionExtraStiffness) + (boneMap(x, modelDepth - y) * boneExtraStiffness)');
	else
		model.physics('solid').feature('lemm1').set('E', 1, 'basalStiffness+(stiffnessMap(x,modelDepth-y)*lesionExtraStiffness)');
	end

	model.physics('solid').feature('lemm1').set('nu_mat', 1, 'userdef');
	model.physics('solid').feature('lemm1').set('nu', 1, 'poissonsRatio');
	model.physics('solid').feature('lemm1').set('rho_mat', 1, 'userdef');
	model.physics('solid').feature('lemm1').set('rho', 1, 'density');

	% set the boundary conditions
	model.physics('solid').feature.create('fix1', 'Fixed', 1);
	model.physics('solid').feature('fix1').selection.set([2]);
	model.physics('solid').feature.create('disp1', 'Displacement1', 1);
	model.physics('solid').feature('disp1').selection.set([3]);
	model.physics('solid').feature('disp1').set('Direction', 2, '1');
	model.physics('solid').feature('disp1').set('U0', 2, '-compression');

	% create the mesh
	model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
	model.mesh('mesh1').feature('size').set('hauto', '1');
	model.mesh('mesh1').run;

	% setup the study and run it
	model.sol.create('sol1');
	model.sol('sol1').study('std1');
	model.sol('sol1').feature.create('st1', 'StudyStep');
	model.sol('sol1').feature('st1').set('study', 'std1');
	model.sol('sol1').feature('st1').set('studystep', 'stat');
	model.sol('sol1').feature.create('v1', 'Variables');
	model.sol('sol1').feature('v1').set('control', 'stat');
	model.sol('sol1').feature.create('s1', 'Stationary');
	model.sol('sol1').feature('s1').feature.create('fc1', 'FullyCoupled');
	model.sol('sol1').feature('s1').feature.remove('fcDef');
	model.sol('sol1').attach('std1');

	model.sol('sol1').runAll;
end