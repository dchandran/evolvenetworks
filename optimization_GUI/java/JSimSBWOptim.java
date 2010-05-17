// JSim Optimizer API for use with SBW

import java.io.*;
import JSim.util.*;
import JSim.data.*;
import JSim.aserver.*;

public class JSimSBWOptim implements ASServer.Messenger {

	private String SBML = null;	
	public void loadSBML(String sbmlText)
	{
		SBML = sbmlText;
	}
	
	private double timeEnd = 100.0;
	public void setEndTime(double t)
	{
		timeEnd = t;
	}
	
	private double stepSize = 0.1;
	public void setStepSize(double dt)
	{
		stepSize = dt;
	}
	
	private String[] targetParams = new String[0];
	public void setTargetParameters(String[] paramNames)
	{
		targetParams = paramNames;
	}
	
	private double[] minParamValue = new double[0];
	public void setMinimumParameterValues(double[] values)
	{
		minParamValue = values;
	}
	
	private double[] maxParamValue = new double[0];
	public void setMaximumParameterValues(double[] values)
	{
		maxParamValue = values;
	}
	
	private String[] targetVars = new String[0];
	public void setTargetVariables(String[] varNames)
	{
		targetVars = varNames;
	}
	
	private String algorithmName = "ggopt";
	public void setAlgorithm(String name)
	{
		algorithmName = name;
	}
	
	private String modelName = "temp";
	public void setModelName(String name)
	{
		modelName = name;
	}
	
	private String targetData = "";
	public void setTargetData(String filename)
	{
		targetData = filename;
	}
	
	private String[] pointWgts = new String[0];
	public void setPointWeights(String[] weights)
	{
		pointWgts = weights;
	}
	
	private int maxCalls = 50;
	public void setMaxCalls(int n)
	{
		maxCalls = n;
	}
	
	private double errTol = 1e-6;
	public void setErrorTolerance(double e)
	{
		errTol = e;
	}
	
	private double stepTol = 1e-6;
	public void setStepTolerance(double e)
	{
		stepTol = e;
	}
	
	private boolean calcCovMat = false;
	public void calculateCovarianceMatrix(boolean b)
	{
		calcCovMat = b;
	}
	
	private boolean savLog = true;
	public void saveLogs(boolean b)
	{
		savLog = b;
	}
	
	private double gradTol = 1e-4;
	public void setGradTolerance(double e)
	{
		gradTol = e;
	}
	
	private double eps = 1e-7;
	public void setEpsilon(double e)
	{
		eps = e;
	}
	
	private int maxIters = 50;
	public void setMaxIterations(int n)
	{
		maxIters = n;
	}
	
	private OptimResults results = null;
	public double[] run() throws Exception
	{
		if (targetParams.length < 1 || 
			targetVars.length < 1)
			return new double[0];
		
		// process command line
	    File sbmlFile = new File(SBML);
	    File csvFile = new File(targetData);
		
		if (!sbmlFile.canRead() || !csvFile.canRead()) return new double[0];
		
	    // initialize connection to gildor (NSR backup server)
	    NamedVal.NList sopts = new NamedVal.NList();
	    sopts.add(NamedVal.create("server", "gamgee.bioeng.washington.edu"));
	    ASServer server = ASServer.create(sopts, this, null);
	    
	    // create model, compile SBML code
	    ASModel model = server.newModelRT();
	    ASInfo.Build build = new ASInfo.Build();
	    build.name = modelName;
	    String sbmlText = UtilIO.readText(sbmlFile);
	    build.modelSource = server.createMML(sbmlText, null);
	    build.sourceType = ASModel.SRC_MML;
	    build.options = new NamedVal.NList();
	    model.buildRT(build);

	    // read reference data file
	    FileReader rdr = new FileReader(csvFile);
	    CSVDataReader crdr = new CSVDataReader(new CSVDataFormat(), rdr);
	    Data.List dataList = crdr.readData();
	    Data refData = dataList.data(0);

	    // create info for Optimizer 
	    ASInfo.Optim oinfo = new ASInfo.Optim(targetParams.length, targetVars.length); // 1 parm, 1 curve to match
	    oinfo.args.alg = algorithmName;
	    
		for (int i=0; i < targetParams.length; ++i)
		{
			String parmName = targetParams[i];
			ASVar parm = model.getASVar(parmName);
			double parmStart = Util.toDouble(parm.getDefault());
			oinfo.args.xname[i] = parmName;
			oinfo.args.xstart[i] = parmStart;
			
			if (minParamValue.length > i)
				oinfo.args.xmin[i] = minParamValue[i];
			else
				oinfo.args.xmin[0] = Double.NaN;
			
			if (maxParamValue.length > i)
				oinfo.args.xmax[i] = maxParamValue[i];
			else
				oinfo.args.xmax[0] = Double.NaN;
		}
				
		oinfo.args.maxCalls = maxCalls;
		oinfo.args.errTol = errTol;
		oinfo.args.stepTol = stepTol;
		oinfo.args.saveLogs = savLog;
	    oinfo.args.calcCovMat = calcCovMat;
	    oinfo.args.maxIters = maxIters;
	    oinfo.args.gradTol = gradTol;
	    oinfo.args.eps = eps;
	    oinfo.baseVals = server.getSolverDefaults();
	    oinfo.refData[0] = refData;
		
		for (int i=0; i < targetVars.length; ++i)
		{
			oinfo.matchExprs[i] = targetVars[i];
			
			if (pointWgts.length > i)
				oinfo.pointWgts[i] = pointWgts[i];
			else
				oinfo.pointWgts[i] = "1";
			oinfo.curveWgts[i] = 1;
		}
		
		String tmax = String.valueOf(timeEnd);
	    String tdelta = String.valueOf(stepSize);
		
	    // run optimizer
	    model.getASVar("time.max").setAssign(tmax);
	    model.getASVar("time.delta").setAssign(tdelta);
	    model.optimRun(oinfo);
	    
	    // get/print optimization results
	    results = model.optimResults();
	    //System.out.println("Final " + parmName + "=" + results.bestX[0]);
	    return results.bestX;
	}

	// ASServer.Messenger implementation (message f. server)
	public void message(ASInfo.Message msg) 
	{
	    System.err.println(msg.text);
	}
}
