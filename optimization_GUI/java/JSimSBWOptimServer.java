import edu.caltech.sbw.*;               // Import SBW.

public class JSimSBWOptimServer
{
    public static void main(String[] args)
    {
        try
        {
            ModuleImpl moduleImpl = new ModuleImpl("JSim Optimization Algorithms");

            moduleImpl.addService("Optimization", "JSim Optimization Algorithms","JSim", JSimSBWOptim.class);
            moduleImpl.run(args);
        }
        catch (SBWException e)
        {
            e.handleWithDialog();
        }
    }
}
