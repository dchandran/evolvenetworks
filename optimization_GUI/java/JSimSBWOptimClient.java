
import edu.caltech.sbw.*;               // Import SBW.

public class JSimSBWOptimClient
{
    public static double[] optimize(String[] args)
    {
		if (args.length < 4)
			return new double[0];

        try
        {
            SBW.connect();

            Module module = SBW.getModuleInstance("JSimSBWOptimServer");
            Service service = module.findServiceByName("Optimization");
            JSimSBWOptimInterface jsimOptim = (JSimSBWOptimInterface) service.getServiceObject(JSimSBWOptimInterface.class);

			jsimOptim.loadSBML(args[0]);
			jsimOptim.setTargetData(args[1]);
			jsimOptim.setModelName(args[2]);
			jsimOptim.setAlgorithm(args[3]);
			jsimOptim.setEndTime(Double.parseDouble(args[4]));
			jsimOptim.setStepSize(Double.parseDouble(args[5]));
			jsimOptim.setErrorTolerance(Double.parseDouble(args[6]));
			jsimOptim.setStepTolerance(Double.parseDouble(args[7]));
			jsimOptim.setGradTolerance(Double.parseDouble(args[8]));
			jsimOptim.setEpsilon(Double.parseDouble(args[9]));
			jsimOptim.setMaxIterations(Integer.parseInt(args[10]));
			
			//jsimOptim.calculateCovarianceMatrix(boolean b);
			//jsimOptim.saveLogs(boolean b);
			
			int n = 0, j = 11;
			for (j=11; j < args.length; ++j, ++n) 
				if (args[j].equals(new String("0")))
					break;
			
			String[] varNames = new String[n];
			j = 11;
			for (int i=0; i < n; ++i, ++j) 
				varNames[i] = args[j];
			++j;
			
			n = (args.length - j)/3;
			String[] params = new String[n];
			double[] minParam = new double[n];
			double[] maxParam = new double[n];
			
			for (int i=0; i < n; ++i, j+=3)
			{
				params[i] = args[j];
				minParam[i] = Double.parseDouble(args[j+1]);
				maxParam[i] = Double.parseDouble(args[j+2]);
			}
			
			jsimOptim.setTargetParameters(params);
			jsimOptim.setTargetVariables(varNames);
			jsimOptim.setMinimumParameterValues(minParam);
			jsimOptim.setMaximumParameterValues(maxParam);
			
			double[] result = jsimOptim.run();

            module.shutdown();
            SBW.disconnect();

            return result;
        }
        catch (SBWException e)
        {
            e.handleWithDialog();
        }

        return new double[0];
    }

    public static void main(String[] args)
    {        
		//"test1.sbml","ref.csv","beta","PX"
		double[] result = optimize(args);
		
		for (int i=0; i < result.length; ++i)
			System.out.print( String.valueOf(result[i]) + " ");

		System.out.println();
        System.exit(0);
    }
}
