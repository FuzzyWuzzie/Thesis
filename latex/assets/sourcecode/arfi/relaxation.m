function model = relaxation(parameters)
	% import COMSOL features
	import com.comsol.model.*
	import com.comsol.model.util.*

	% create the model
	model = ModelUtil.create('Model');
	model.modelPath('ARFI');
	model.name('relaxation.mph');

	% set parameters in the model
	model.param.set('domainWidth', sprintf('%f[m]', parameters.domainWidth));
	model.param.set('domainDepth', sprintf('%f[m]', parameters.domainDepth));
	model.param.set('focalDepth', sprintf('%f[mm]', parameters.focalDepth));
	model.param.set('timeStep', sprintf('%f[us]', parameters.timeStep));
	model.param.set('loadTime', sprintf('%f[us]', parameters.loadTime));
	model.param.set('listenTime', sprintf('%f[ms]', parameters.listenTime));
	model.param.set('focalY', 'domainDepth-focalDepth');
	model.param.set('dFocalY', 'focalY/1[m]');
	model.param.set('stiffnessRatio', sprintf('%f', parameters.stiffnessRatio));
	model.param.set('cutoffAmplitude', sprintf('%f', parameters.cutoffAmplitude));

	% load interpolation files for both geometry
	% and initial acoustic radiation force
	model.modelNode.create('mod1');
	model.func.create('int1', 'Interpolation');
	model.func.create('int2', 'Interpolation');
	model.func.create('step1', 'Step');
	model.func.create('im1', 'Image');
	model.func('int1').set('funcs', {'Fx' '1'});
	model.func('int1').set('source', 'file');
	model.func('int1').set('filename', 'ARFI/Fx.txt');
	model.func('int1').set('struct', 'grid');
	model.func('int1').set('defvars', 'on');
	model.func('int1').set('extrap', 'value');
	model.func('int1').set('argunit', 'm');
	model.func('int1').set('fununit', 'N/(m^3)');
	model.func('int2').set('funcs', {'Fy' '1'});
	model.func('int2').set('source', 'file');
	model.func('int2').set('filename', 'ARFI/Fy.txt');
	model.func('int2').set('struct', 'grid');
	model.func('int2').set('defvars', 'on');
	model.func('int2').set('extrap', 'value');
	model.func('int2').set('argunit', 'm');
	model.func('int2').set('fununit', 'N/(m^3)');
	model.func('step1').set('location', 'loadTime');
	model.func('step1').set('from', '1');
	model.func('step1').set('to', '0');
	model.func('step1').set('smoothactive', false);
	model.func('im1').set('funcname', 'stiffnessMap');
	model.func('im1').set('filename', 'ARFI/stiffnessMap.png');
	model.func('im1').set('xmin', '-domainWidth/2');
	model.func('im1').set('xmax', 'domainWidth/2');
	model.func('im1').set('ymax', 'domainDepth');
	model.func('im1').set('extrap', 'value');

	% create the geometry
	model.geom.create('geom1', 2);
	model.geom('geom1').feature.create('r1', 'Rectangle');
	model.geom('geom1').feature.create('r2', 'Rectangle');
	model.geom('geom1').feature('r1').set('pos', {'0' '0'});
	model.geom('geom1').feature('r1').set('size', {'domainWidth/2' 'domainDepth'});
	model.geom('geom1').feature('r2').set('pos', {'domainWidth/2' '0'});
	model.geom('geom1').feature('r2').set('size', {'domainWidth/8' 'domainDepth'});
	model.geom('geom1').run;

	% set up material properties
	model.material.create('mat1');
	model.material('mat1').propertyGroup.create('KG', 'Bulk modulus and shear modulus');
	model.material('mat1').propertyGroup('def').set('density', '1060');
	model.material('mat1').propertyGroup('KG').set('K', '');
	model.material('mat1').propertyGroup('KG').set('G', '');
	model.material('mat1').propertyGroup('KG').set('K', 'stiffnessRatio*515.656[kPa]');
	model.material('mat1').propertyGroup('KG').set('G', 'stiffnessRatio*1032[Pa]');

	% set up the physics
	% including boundary and ``initial'' conditions
	% fix the bottom of the domain
	model.physics.create('solid', 'SolidMechanics', 'geom1');
	model.physics('solid').feature.create('fix1', 'Fixed', 1);
	model.physics('solid').feature('fix1').selection.set([2 5]);
	model.physics('solid').feature('fix1').name('bottom hold');

	% prevent motion in the vertical direction at the top boundary
	model.physics('solid').feature.create('disp1', 'Displacement1', 1);
	model.physics('solid').feature('disp1').selection.set([3 6]);
	model.physics('solid').feature('disp1').set('Direction', {'0'; '1'; '0'});
	model.physics('solid').feature('disp1').name('top hold');

	% the ARFI is a body load
	model.physics('solid').feature.create('bl1', 'BodyLoad', 2);
	model.physics('solid').feature('bl1').selection.set([1]);
	model.physics('solid').feature('bl1').set('FperVol', {'Fx(x,y)*step1(t[1/s])'; 'Fy(x,y)*step1(t[1/s])'; '0'});
	model.physics('solid').feature('bl1').set('Ftot', {''; 'Fy(x,y)*step1(t[1/s])*10^2'; '0'});
	model.physics('solid').feature('bl1').name('ARFI load');

	% make the model symmetric
	model.physics('solid').feature.create('sym1', 'SymmetrySolid', 1);
	model.physics('solid').feature('sym1').selection.set([1]);

	% use a viscoelastic tissue model
	model.physics('solid').feature.create('vmm1', 'ViscoelasticModel', 2);
	model.physics('solid').feature('vmm1').selection.all;
	model.physics('solid').feature('vmm1').set('K_mat', 'userdef');
	model.physics('solid').feature('vmm1').set('K', '(515.656 [kPa]) * (1 + (stiffnessMap(x, y) * (stiffnessRatio - 1)))');
	model.physics('solid').feature('vmm1').set('G_mat', 'userdef');
	model.physics('solid').feature('vmm1').set('G', '(1032 [Pa]) * (1 + (stiffnessMap(x, y) * (stiffnessRatio - 1)))');
	model.physics('solid').feature('vmm1').set('Branch', {'1'; '2'; '2'});
	model.physics('solid').feature('vmm1').set('Gi', {'(791 [Pa]) * (1 + (stiffnessMap(x, y) * (stiffnessRatio - 1)))'; '(66.5 [Pa]) * (1 + (stiffnessMap(x, y) * (stiffnessRatio - 1)))'; '(628 [mPa]) * (1 + (stiffnessMap(x, y) * (stiffnessRatio - 1)))'});
	model.physics('solid').feature('vmm1').set('tau', {'2 [s]'; '40 [s]'; '80 [s]'});
	model.physics('solid').feature('vmm1').set('NearlyIncompressible', '1');

	% create the mesh
	model.mesh.create('mesh1', 'geom1');
	model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
	model.mesh('mesh1').feature.create('map1', 'Map');
	model.mesh('mesh1').feature('size').set('hauto', 3);
	model.mesh('mesh1').feature('map1').active(false);
	model.mesh('mesh1').feature('map1').set('adjustedgdistr', true);
	model.mesh('mesh1').run;

	% set up a probe at the focal point
	model.result.table.create('tbl1', 'Table');
	model.result.table.create('tbl2', 'Table');
	model.probe.create('pdom1', 'DomainPoint');
	model.probe('pdom1').model('mod1');
	model.coordSystem.create('pml1', 'geom1', 'PML');
	model.coordSystem('pml1').selection.set([2]);
	model.probe('pdom1').set('coords2', {'0' 'domainDepth-focalDepth'});
	model.probe('pdom1').feature('ppb1').set('probename', 'focalPointDisplacement');
	model.probe('pdom1').feature('ppb1').set('table', 'tbl1');
	model.probe('pdom1').feature('ppb1').set('window', 'window1');
	model.result.table('tbl1').name('Probe Table 1');
	model.result.table('tbl2').comments('Global Evaluation 1 (t)');

	% create a transient study
	model.study.create('std1');
	model.study('std1').feature.create('time', 'Transient');

	% set up the solution parameters
	model.sol.create('sol1');
	model.sol('sol1').study('std1');
	model.sol('sol1').attach('std1');
	model.sol('sol1').feature.create('st1', 'StudyStep');
	model.sol('sol1').feature.create('v1', 'Variables');
	model.sol('sol1').feature.create('t1', 'Time');
	model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
	model.sol('sol1').feature('t1').feature.create('st1', 'StopCondition');
	model.sol('sol1').feature('t1').feature.remove('fcDef');
	model.sol('sol1').attach('std1');
	model.sol('sol1').feature('st1').name('Compile Equations: Time Dependent');
	model.sol('sol1').feature('st1').set('studystep', 'time');
	model.sol('sol1').feature('v1').set('control', 'time');
	model.sol('sol1').feature('v1').feature('mod1_u').set('scalemethod', 'manual');
	model.sol('sol1').feature('v1').feature('mod1_u').set('scaleval', '1e-2*0.0640312423743285');
	model.sol('sol1').feature('t1').set('tlist', 'range(0,timeStep,loadTime+listenTime)');
	model.sol('sol1').feature('t1').set('fieldselection', 'mod1_u');
	model.sol('sol1').feature('t1').set('atolmethod', {'mod1_u' 'global' 'mod1_solid_qXX3' 'global' 'mod1_solid_qXY2' 'global' 'mod1_solid_qXY3' 'global' 'mod1_solid_qXX1' 'global'  ...
	'mod1_solid_qXX2' 'global' 'mod1_solid_qYY3' 'global' 'mod1_solid_qYY2' 'global' 'mod1_solid_qYY1' 'global' 'mod1_solid_qXY1' 'global'  ...
	'mod1_solid_pw' 'global'});
	model.sol('sol1').feature('t1').set('atol', {'mod1_u' '1e-3' 'mod1_solid_qXX3' '1e-3' 'mod1_solid_qXY2' '1e-3' 'mod1_solid_qXY3' '1e-3' 'mod1_solid_qXX1' '1e-3'  ...
	'mod1_solid_qXX2' '1e-3' 'mod1_solid_qYY3' '1e-3' 'mod1_solid_qYY2' '1e-3' 'mod1_solid_qYY1' '1e-3' 'mod1_solid_qXY1' '1e-3'  ...
	'mod1_solid_pw' '1e-3'});
	model.sol('sol1').feature('t1').set('atoludot', {'mod1_u' '1e-3' 'mod1_solid_qXX3' '1e-3' 'mod1_solid_qXY2' '1e-3' 'mod1_solid_qXY3' '1e-3' 'mod1_solid_qXX1' '1e-3'  ...
	'mod1_solid_qXX2' '1e-3' 'mod1_solid_qYY3' '1e-3' 'mod1_solid_qYY2' '1e-3' 'mod1_solid_qYY1' '1e-3' 'mod1_solid_qXY1' '1e-3'  ...
	'mod1_solid_pw' '1e-3'});
	model.sol('sol1').feature('t1').set('atoludotactive', {'mod1_u' 'off' 'mod1_solid_qXX3' 'off' 'mod1_solid_qXY2' 'off' 'mod1_solid_qXY3' 'off' 'mod1_solid_qXX1' 'off'  ...
	'mod1_solid_qXX2' 'off' 'mod1_solid_qYY3' 'off' 'mod1_solid_qYY2' 'off' 'mod1_solid_qYY1' 'off' 'mod1_solid_qXY1' 'off'  ...
	'mod1_solid_pw' 'off'});
	model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
	model.sol('sol1').feature('t1').set('tstepsgenalpha', 'manual');
	model.sol('sol1').feature('t1').set('timestepgenalpha', 'timeStep');
	model.sol('sol1').feature('t1').set('plot', 'on');
	model.sol('sol1').feature('t1').set('plotgroup', 'pg2');
	model.sol('sol1').feature('t1').feature('fc1').set('plot', 'on');
	model.sol('sol1').feature('t1').feature('fc1').set('plotgroup', 'pg2');

	% add a stop condition to stop the simulation when the tissue has relaxed
	model.sol('sol1').feature('t1').feature('st1').set('stopcond', 'if(t > 0.003 && mod1.focalPointDisplacement < cutoffAmplitude, -1, 1)');

	% create data sets to evaluate later:
	%	* The displacement at the focal point over time
	%	* The displacement along an axial line going through the focal point
	%	* The displacement along a lateral line going through the focal point
	model.result.dataset.create('cpt1', 'CutPoint2D');
	model.result.dataset.create('cln1', 'CutLine2D');
	model.result.dataset.create('cln2', 'CutLine2D');
	model.result.dataset.create('dset2', 'Solution');
	model.result.dataset('dset2').set('probetag', 'pdom1');
	model.result.dataset.create('cpt2', 'CutPoint2D');
	model.result.dataset('cpt2').set('probetag', 'pdom1');
	model.result.dataset('cpt2').set('data', 'dset2');
	model.result.numerical.create('pev1', 'EvalPoint');
	model.result.numerical('pev1').set('probetag', 'ppb1');
	model.result.numerical.create('gev1', 'EvalGlobal');
	model.result.numerical('gev1').set('probetag', 'none');
	model.result.create('pg1', 'PlotGroup1D');
	model.result('pg1').set('probetag', 'none');
	model.result('pg1').feature.create('ptgr1', 'PointGraph');
	model.result.create('pg2', 'PlotGroup2D');
	model.result('pg2').feature.create('surf1', 'Surface');
	model.result.create('pg3', 'PlotGroup1D');
	model.result('pg3').set('probetag', 'none');
	model.result('pg3').feature.create('lngr1', 'LineGraph');
	model.result.create('pg4', 'PlotGroup1D');
	model.result('pg4').set('probetag', 'none');
	model.result('pg4').feature.create('lngr1', 'LineGraph');
	model.result.create('pg5', 'PlotGroup1D');
	model.result('pg5').set('probetag', 'window1');
	model.result('pg5').feature.create('tblp1', 'Table');
	model.result('pg5').feature('tblp1').set('probetag', 'ppb1');
	model.result.dataset('cpt1').name('Focal Point');
	model.result.dataset('cpt1').set('pointx', '0');
	model.result.dataset('cpt1').set('pointy', 'dFocalY');
	model.result.dataset('cln1').name('Lateral Cut');
	model.result.dataset('cln1').set('genpoints', {'-domainWidth/2' 'domainDepth - focalDepth'; 'domainWidth/2' 'domainDepth - focalDepth'});
	model.result.dataset('cln2').name('Axial Cut');
	model.result.dataset('cln2').set('genpoints', {'0' 'domainDepth'; '0' '0'});
	model.result.dataset('dset2').name('Probe Solution 2');
	model.result.dataset('cpt2').set('pointy', 'focalY');
	model.result.numerical('gev1').set('table', 'tbl2');
	model.result.numerical('gev1').set('expr', 't');
	model.result.numerical('gev1').set('unit', 's');
	model.result.numerical('gev1').set('descr', 'Time');
	model.result.numerical('gev1').set('dataseries', 'maximum');
	model.result.numerical('pev1').setResult;
	model.result.numerical('gev1').setResult;
	model.result('pg1').name('Focal Point Relaxation');
	model.result('pg1').set('data', 'cpt1');
	model.result('pg1').set('xlabel', 'Time (ms)');
	model.result('pg1').set('ylabel', 'Displacement (m)');
	model.result('pg1').set('xlabelactive', false);
	model.result('pg1').set('ylabelactive', false);
	model.result('pg1').feature('ptgr1').set('descractive', true);
	model.result('pg1').feature('ptgr1').set('descr', 'Displacement');
	model.result('pg1').feature('ptgr1').set('titletype', 'manual');
	model.result('pg1').feature('ptgr1').set('title', 'Focal Point Relaxation');
	model.result('pg1').feature('ptgr1').set('xdata', 'expr');
	model.result('pg1').feature('ptgr1').set('xdataexpr', 't');
	model.result('pg1').feature('ptgr1').set('xdataunit', 'ms');
	model.result('pg1').feature('ptgr1').set('xdatadescr', 'Time');
	model.result('pg2').name('Surface Displacement');
	model.result('pg3').name('Displacement of Axial Focal Cut');
	model.result('pg3').set('data', 'cln2');
	model.result('pg3').set('xlabel', 'Depth (m) (m)');
	model.result('pg3').set('ylabel', 'Total displacement (m)');
	model.result('pg3').set('xlabelactive', false);
	model.result('pg3').set('ylabelactive', false);
	model.result('pg3').feature('lngr1').set('xdata', 'expr');
	model.result('pg3').feature('lngr1').set('xdataexpr', 'domainDepth - y');
	model.result('pg3').feature('lngr1').set('xdatadescractive', true);
	model.result('pg3').feature('lngr1').set('xdatadescr', 'Depth (m)');
	model.result('pg3').feature('lngr1').set('legend', true);
	model.result('pg4').name('Displacement of Lateral Focal Cut');
	model.result('pg4').set('data', 'cln1');
	model.result('pg4').set('xlabel', 'x-coordinate (m)');
	model.result('pg4').set('ylabel', 'abs(v) (m)');
	model.result('pg4').set('legendpos', 'lowerright');
	model.result('pg4').set('xlabelactive', false);
	model.result('pg4').set('ylabelactive', false);
	model.result('pg4').feature('lngr1').set('expr', 'abs(v)');
	model.result('pg4').feature('lngr1').set('descr', 'abs(v)');
	model.result('pg4').feature('lngr1').set('xdata', 'expr');
	model.result('pg4').feature('lngr1').set('xdataexpr', 'x');
	model.result('pg4').feature('lngr1').set('xdatadescr', 'x-coordinate');
	model.result('pg4').feature('lngr1').set('legend', true);
	model.result('pg5').name('Probe 1D Plot Group 5');
	model.result('pg5').set('xlabel', 't');
	model.result('pg5').set('ylabel', 'Total displacement, Point Probe Expression 1');
	model.result('pg5').set('windowtitle', 'Probe Plot 1');
	model.result('pg5').set('xlabelactive', false);
	model.result('pg5').set('ylabelactive', false);
	model.result('pg5').feature('tblp1').name('Probe Table Graph 1');

	% initiate the domain probe we defined earlier
	model.probe('pdom1').genResult([]);

	% set the timestepping for the solution
	model.study('std1').feature('time').set('tlist', 'range(0,timeStep,loadTime+listenTime)');
	model.study('std1').feature('time').set('plot', 'on');
	model.study('std1').feature('time').set('plotgroup', 'pg2');

	% run the model
	model.sol('sol1').runAll;
end
