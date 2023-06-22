import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.LinkedHashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Comparator;
import java.util.Collections;

class Vertex {

  float x;
  float y;
  float z;
  
  int pathNo;
   
  public Vertex() {}

  public Vertex(float x, float y, float z) {

    this.x = x;
    this.y = y;
    this.z = z;
    this.pathNo = 999;
  }
  public Vertex(double x, double y, double z) {

    this.x = (float)x;
    this.y = (float)y;
    this.z = (float)z;
    this.pathNo = 999;
  }
  public static Vertex sub(Vertex v1,Vertex v2) {
    float x = v1.x - v2.x;
    float y = v1.y - v2.y;
    float z = v1.z - v2.z;
    return new Vertex(x,y,z);

  }
  public static Vertex add(Vertex v1, Vertex v2) {
    return new Vertex(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z); 
  }
  public static Vertex add(Vertex v1, Vertex v2, float k) {
    return new Vertex(v1.x + v2.x*k, v1.y + v2.y*k, v1.z + v2.z*k); 
  }
  public String toString() {
    return "( " + x + "," + y + "," + z + " )";
  }

  //@Override
  //public boolean equals(Object o) {
    
  //}
}
class Line {
  Vertex v1;
  Vertex v2;
  
  public Line(Vertex v1, Vertex v2) {

    this.v1 = v1;
    this.v2 = v2;

  }
  public Line(Line line) {
    this.v1 = line.v1;
    this.v2 = line.v2;
  }
  public Vertex getVertex1() {
    return v1;
  }
  public Vertex getVertex2() {
    return v2;
  }
  public Line reverseLine() {
    return new Line(v2, v1);
  }
  public void swapVerts() {
    Vertex t = v1;
    v1 = v2;
    v2 = t;
  }
  public double distanceVertex1(Line l) {
    double dx = l.v1.x - v1.x;
    double dy = l.v1.y - v1.y;
    double d = dx * dx + dy * dy;
    return d;
  }
  public double distanceVertex2(Line l) {
    double dx = l.v2.x - v1.x;
    double dy = l.v2.y - v1.y;
    double d = dx * dx + dy * dy;
    return d;
  }
  public String toString() {
    return "Line: " + v1 + " to " + v2;
  }
  public static Line add(Line a, Line b) {
    
    Vertex p1 = Vertex.add(a.getVertex1(),b.getVertex1());
    Vertex p2 = Vertex.add(a.getVertex2(),b.getVertex2());
    return new Line(p1,p2);

  }
  //unit vector in Vertex b
  public static Line add(Line a, Vertex b) {
    
    Vertex p1 = Vertex.add(a.getVertex1(),b);
    Vertex p2 = Vertex.add(a.getVertex2(),b);
    return new Line(p1,p2);

  }

  public static Line add(Line a, Vertex b, float k) {
    
    Vertex p1 = Vertex.add(a.getVertex1(),b,k);
    Vertex p2 = Vertex.add(a.getVertex2(),b,k);
    return new Line(p1,p2);

  }

  //unit vector in Vertex b
  public Vertex getNormal(){

            double x0 = getVertex1().x;
            double x1 = getVertex2().x;

            double y0 = getVertex1().y;
            double y1 = getVertex2().y;

            double ux = x1-x0;
            double uy = y1-y0;
            double rx = 0;
            double ry = 1;
            double nx = ux*rx - uy*ry;
            double ny = ux*ry + uy*rx;
            nx = -uy;
            ny =  ux;
            double mag = (double)Math.sqrt(nx*nx + ny*ny);
            nx /= mag;
            ny /= mag;
            //System.out.println("normal: " + nx + " " + ny);
            return new Vertex((float)nx,(float)ny,0);
            
  }

  public static Line mergeLines(Line a, Line b) {
    Vertex v1 = a.getVertex1();
    Vertex v2 = b.getVertex2();
    return new Line(v1,v2);
  }
  public static boolean isColinear(Line a, Line b) {
      float x1 = a.getVertex1().x;
            float x2 = a.getVertex2().x;

            float y1 = a.getVertex1().y;
            float y2 = a.getVertex2().y;

            float x3 = b.getVertex1().x;
            float x4 = b.getVertex2().x;

            float y3 = b.getVertex1().y;
            float y4 = b.getVertex2().y;

            double d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

            if(Math.abs(d) > 0.0001) return false;
            else return true;

  }
  public static Vertex getIntersection(Line a,Line b,boolean c){
            
            //System.out.println("(Line A: " + a + ")");
            //System.out.println("(Line B: " + b + ")");

            float x1 = a.getVertex1().x;
            float x2 = a.getVertex2().x;

            float y1 = a.getVertex1().y;
            float y2 = a.getVertex2().y;

            float x3 = b.getVertex1().x;
            float x4 = b.getVertex2().x;

            float y3 = b.getVertex1().y;
            float y4 = b.getVertex2().y;

            double d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
            double px = 0.0;
            double py = 0.0;
            
            //if (d!=0.0) {
            if(Math.abs(d) > 0.0001) {
                 px = (x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4);
                 px = px/d;
                 py = (x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4);
                 py = py/d; 
            } 
            else {
                 //System.out.println("coincident lines: fix this algorithm");

                 px = x1+(x4-x2);  
                 //px = x2;
                 py = y1+(y4-y2);
                 //py = y2;  
                 
            }

            return new Vertex((float)px,(float)py,0);

  }
  public static Vertex getIntersection(Line a,Line b){
            
            float ax0 = a.getVertex1().x;
            float ax1 = a.getVertex2().x;

            float ay0 = a.getVertex1().y;
            float ay1 = a.getVertex2().y;

            float bx0 = b.getVertex1().x;
            float bx1 = b.getVertex2().x;

            float by0 = b.getVertex1().y;
            float by1 = b.getVertex2().y;

            float dxa = ax1-ax0;
            float dxb = bx1-bx0;
            System.out.println(a + " " + b);

            if(dxa == 0.0f && dxb==0.0f){
                System.out.println("2 infinite slope");
                System.out.println(a + " " + b);
                float x = ax0 + 0.5f*(bx1-ax0);
                float y = ay0 + 0.5f*(by1-ay0);
                return new Vertex(x,y,0);
            }
            if(dxa == 0.0f){
                
                float mb = (by1-by0)/dxb;
                float binter = by1 - mb*bx1;
                float x = ax0;
                float y = by0 + mb*(by1-by0);
                System.out.println("a slope inf: " + x + " " + y);
                return new Vertex(x,y,0);
            }
            if(dxb == 0.0f){
                
                float ma = (ay1-ay0)/dxa;
                float ainter = ay1 - ma*ax1;
                float x = bx0 ;
                float y = ay0 + ma*(ax1-ax0);
                System.out.println("b slope inf" + x + " " + y);
                return new Vertex(x,y,0);
            }

            //clip infinitely long lines
            double ma = (ay1-ay0)/(ax1-ax0);
            double mb = (by1-by0)/(bx1-bx0);
            double ainter = ay1-ma*ax1;
            double binter = by1-mb*bx1;
            
            System.out.println("a param: " + ma + "  " + ainter + (ma<-Float.MAX_VALUE || ma>Float.MAX_VALUE));
            System.out.println("b param: " + mb + "  " + binter + (mb<-Float.MAX_VALUE || mb>Float.MAX_VALUE));

          
            // calc intercept by substitution
            double d = (binter-ainter)/(ma-mb);
            double f = (binter * ma - ainter*ma)/(ma-mb) + ainter;
            double y = f;
            double x = d;

            double dm = ma-mb;
            dm = Math.abs(dm);

            if(ma==mb || dm < 1e-3) {
              System.out.println("slope equal");
              x = ax0 + 0.5f*(bx1-ax0);
              y = ay0 + 0.5f*(by1-ay0);
              return new Vertex(x,y,0); 
            }
            if(ma==0.0f) {
             y = ainter;
             x = (ainter-binter)/mb;
              
            }
            System.out.println("x,y "+ x + " "+ y);
            return new Vertex(x,y,0);

  }
  public Vertex findIntersection(Line b){
            float x1 = this.v1.x;
            float x2 = this.v2.x;
            float y1 = this.v1.y;
            float y2 = this.v2.y;
            float x3 = b.v1.x;
            float x4 = b.v2.x;
            float y3 = b.v1.y;
            float y4 = b.v2.y;

            float t = (x1-x3)*(y3-y4)-(y1-y3)*(x3-x4);
            float td= (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
            
            float u = (x1-x3)*(y1-y2)-(y1-y3)*(x1-x2);
            t= t/td;
            u= u/td;
            
            if( (t>=0 && t<=1) && (u>=0 && u<=1) ) {
              //System.out.println((t));
              //System.out.println((u));
              float px = x1 + t*(x2-x1);
              float py = y1 + t*(y2-y1);
              //System.out.println(px+"  "+py);
              return new Vertex(px,py,0);
            }
            
            return null;  

  }
}

class Tri {

  Vertex v1;
  Vertex v2;
  Vertex v3;

  boolean intersectFound;

  public Tri(Vertex v1, Vertex v2, Vertex v3) {
    this.v1 = v1;
    this.v2 = v2;
    this.v3 = v3;
  }
  public boolean isBetween(float z,float z1,float z2) {
    float temp;
    if(z1 > z2) {
      temp = z1;
      z1 = z2;
      z2 = temp;
    } 
    if(z >= z1 && z <= z2) return true;

    return false;

  }
  public Vertex getNormal() {
    Vertex n = new Vertex();
    Vertex u = Vertex.sub(v2,v1);
    Vertex w = Vertex.sub(v3,v1);

    n.x = u.y * w.z - u.z * w.y;
    n.y = u.z * w.x - u.x * w.z;
    n.z = u.x * w.y - u.y * w.x;

    return n; 

  }
  public boolean zIntersects(float z) {

    if( isBetween(z, v1.z, v2.z) ) return true;
    if( isBetween(z, v2.z, v3.z) ) return true;
    if( isBetween(z, v3.z, v1.z) ) return true;
    return false;

  }

  public static Vertex getInter(float z, Vertex v1, Vertex v2) {
    double mxz;
    
    double myz;
    double intercept;
  
    double x;
    double y;
   
    mxz = (v2.x - v1.x)/(v2.z - v1.z);
    intercept = v1.x - v1.z*mxz;
   
    x = mxz*(z) + intercept;
    
    myz = (v2.y - v1.y)/(v2.z - v1.z);
    intercept = v1.y - v1.z*myz;    

    y = myz*(z) + intercept;
    //System.out.println(v1);
    //System.out.println(v2);
    //System.out.println("mxz: " + mxz);
    //System.out.println("myz: " + myz);
    
    return new Vertex((float)x,(float)y,(float)z);
   
  }

  public Line getIntersects(float z) {

    int i = 0;
    intersectFound = false;
    Vertex[] verts = new Vertex[2];

    //if(v1.z == v2.z && v1.z == v2.z) {
    //    return new Line(null,null);
    //}
    if( isBetween(z, v1.z, v2.z) ) {

      Vertex v = getInter(z,v1,v2);
      
      verts[i] = v;
      ++i;
      intersectFound = true; 
      //System.out.println("E12: " + v);

    }
    if( isBetween(z, v2.z, v3.z) ) {

      Vertex v = getInter(z,v2,v3);
      verts[i] = v;
      ++i;
      intersectFound = true; 
      //System.out.println("E23: " + v);

    }
    if( isBetween(z, v3.z, v1.z) ) {

      Vertex v = getInter(z,v3,v1);
      verts[i] = v;
      ++i; 
      intersectFound = true; 
      //System.out.println("E31: " + v);

    }

    /*
    if(i==0) {
      System.out.println("No intersect");
      System.out.println(z);
      System.out.println(v1);
      System.out.println(v2);
      System.out.println(v3);
    }
    */

    if(i==3) System.out.println("(3 intersects found)");


    return new Line(verts[0],verts[1]);

  }

}
class DistanceLine implements Comparator<DistanceLine> {
  static double distanceMaxError = 0.001; 
  Line line;
  double distance;
  double reverseDistance;
  double dx;
  double dy;

  double rdx;
  double rdy;

  public DistanceLine(Line line, Line compare) {

    this.line = line;

    dx = compare.getVertex2().x - line.getVertex1().x;
    dy = compare.getVertex2().y - line.getVertex1().y;
    
    rdx = compare.getVertex2().x - line.getVertex2().x;
    rdy = compare.getVertex2().y - line.getVertex2().y;
  }
  public DistanceLine(Line line, double distance, double reversedDistance) {
    this.line = line;
    this.distance = distance;
    this.reverseDistance = reverseDistance;
    
  }

  public Line getLine() {
    return line;
  }
  public boolean endsMatch() {
    // absolute value
    if(dx<0) dx=-dx;
    if(dy<0) dy=-dy;
    
    if(dx < distanceMaxError && dy < distanceMaxError) {
      return true;
    }
    return false;
  }

  public boolean endsReverseMatch() {
    // absolute value
    if(rdx<0) rdx=-rdx;
    if(rdy<0) rdy=-rdy;
    
    if(rdx < distanceMaxError && rdy < distanceMaxError) {
      return true;
    }
    return false;
  }
  @Override
  public int compare(DistanceLine l1, DistanceLine l2) {
    double d = l1.distance - l2.distance;
    int r=0;
    if(d==0) r = 0;
    if(d>0) r = 1;
    if(d<0) r = -1;
    return r;
  }
  
  public String toString() {
    String r = line.toString() + "  " + distance + "  " + reverseDistance;
    return r;
  }
}
 
class STLwriter {

  public static void writeFacet( Tri tri ) {

    StringBuilder sb = new StringBuilder();
    sb.append("facet normal ");
    sb.append(tri.getNormal());
    sb.append("\n\touter loop\n");
    sb.append("\t\tvertex " + tri.v1);
    sb.append("\n\t\tvertex " + tri.v2);
    sb.append("\n\t\tvertex " + tri.v3);
    sb.append("\n\tendloop\n");
    sb.append("endfacet\n");

    System.out.println(sb);
  } 
  

}
class SortByY implements Comparator<Vertex> {
  
  public int compare(Vertex a, Vertex b) 
  {
    float r = a.y-b.y;
    if(r==0.0f) return 0;
    if(r < 0.0f) return -1;
    return 1;

  }


}

class ShellBuildDirection {
  public int pathNo;
  public String shellType="";

  public void set(int pathNo, String type) {

    this.pathNo = pathNo;

    if(type.equals("INNER") && shellType.equals("OUTER")) {
        shellType = "INNER";
    }
    else {
        shellType = type;
    }

  }

  public String toString() {
    return pathNo + "  " + shellType;
  }

  @Override
  public boolean equals(Object o) {
    if(o==null) return false;
    if(getClass() != o.getClass()) return false;

    ShellBuildDirection s = (ShellBuildDirection)o;
    if( (pathNo == s.pathNo) && (shellType.equals(s.shellType)) ) return true;

    return false; 

  }

  public int hashCode() {

    return 31 * pathNo + shellType.hashCode();

  }

}

class ShellBuildSet {

    Set<ShellBuildDirection> shellSet = new HashSet<>();

    // An INNER shelltype can override
    // an OUTER shelltype
    // however if an INNER shell type is already set for the path
    // then it can not be set as an OUTER
    public void add(ShellBuildDirection s) {
    
        ShellBuildDirection sbd = new ShellBuildDirection();
        sbd.set(s.pathNo, "OUTER");
       
        if(s.shellType.equals("INNER")) {
          
          if (shellSet.contains(sbd)) {
            shellSet.remove(sbd);
          }
          
          shellSet.add(s);
          
        } else if (s.shellType.equals("OUTER")) {

          sbd.set(s.pathNo, "INNER");
          if(!shellSet.contains(sbd)) {
            shellSet.add(s);
          }

        } 
        else {
          
          shellSet.add(s);
        }
       
    }

    public String findByPathNo(int pathNo){

      for(ShellBuildDirection s : shellSet) {
        if(s.pathNo == pathNo) {
          return s.shellType;
        }
      }
      return "NONE";
    }

    public String toString() {
      return shellSet.toString();
    }

}
public class Triangulate {

    static List<Line> p = new ArrayList<>();

    static ShellBuildSet shell = new ShellBuildSet();

    static float minX = Float.MAX_VALUE;
    static float maxX = -Float.MAX_VALUE;

    static float minY = Float.MAX_VALUE;
    static float maxY = -Float.MAX_VALUE;

    static float minZ = Float.MAX_VALUE;
    static float maxZ = -Float.MAX_VALUE;

    /**
    *
    */
    public static List<Line> findPath(List<Line> lines, int pathIndex, boolean print) {
  
      Line seed = p.get(pathIndex);
      //lines.remove(0);

      Line foundLine = null;
      boolean reverse = false;
      for(Line l : lines) {
        

        DistanceLine dl = new DistanceLine(l,seed);
        if(dl.endsMatch()) {
          //System.out.println("hit: "+l);
          p.add(l);
          foundLine = l;
          break;
        }
        else if (dl.endsReverseMatch()) {
          //System.out.println("rhit: "+l);
          p.add(l.reverseLine());
          foundLine = l;
          reverse = true;
          break;
        }
        
        
      }//for
      
      lines.remove(foundLine);
      //lines.remove(foundLine);
      //lines.remove(foundLine);
      
      if(false) {
        System.out.println("\npath is now: ");
        for(Line l : p) {
          System.out.println(l); 
        }
      }
      
      if(false) {
      
        System.out.println("(lines is now (size = " + lines.size() + ") : )");
        for(Line l : lines) {
          System.out.println("(" + l + ")");
        }
      }

      return lines;
    }
    public static List<Line> buildPath(List<Line> lines, boolean print) {

            //seed the start
            
            p.add(lines.get(0));
            lines.remove(0);
            int s = lines.size();

            for(int j=0;j<s;j++) {

              lines = findPath(lines,p.size()-1,print);

            }
            if(print){
               System.out.println("(path found is (size = " + lines.size() + ") : )");
               for(Line l : p) {
                 System.out.println("(" + l + ")");
               }
            }
            return lines;
    }

    public static List<Line> merge(List<Line> path) {
       List<Line> newPath = new ArrayList<>();
       
       //start
       Line g = path.get(0);
       newPath.add(g);
       int j=0;
       for(int i=1;i<path.size();i++) {
         
         Line lineA = newPath.get(j);
         Line lineB = path.get(i);
         boolean co = Line.isColinear(lineA, lineB);
         if(co) {
           Line m = Line.mergeLines(lineA, lineB);
           newPath.set(j, m);
           
         }else {
           newPath.add(lineB);
           j+=1;
         }
         
           
       }
       System.out.println("( was: "+path.size() + " now: " + j + ")");
       return newPath;
    }
    // Line a and b must to joined in path by a's vertex and b's vertex 2 the same.
    // true == ccw
    public static Vertex getWindiness(Line a, Line b) {
      //System.out.println(a + "  " + b);
      Vertex v1 = a.getVertex1();
      Vertex v2 = a.getVertex2();
      Vertex v3 = b.getVertex2();
      Tri t = new Tri(v1,v2,v3);

      Vertex n = t.getNormal();

      //System.out.println("( windiness normal : " + n + " )");

      //if(n.z>0) return true;
      //else return false;
      return n;
    }
    public static boolean getWindiness(List<Line> path) {
      Line lineA = path.get(0);
      Vertex n = new Vertex(0,0,0);

      for(int i=1;i<path.size();i++) {
        Line lineB = path.get(i);
        Vertex s = getWindiness(lineA, lineB);
        n = Vertex.add(n,s);
        lineA = lineB;
      }
      System.out.println("( windiness normal : " + n + " )");

      if(n.z>0) return true;
      else return false;
      
      
    }
    /** This uses area calc to determine the windiness
    *   positive indicates a clockwise winding
    *   This method returns true indicates counter clockwise
    *   BTW multiply by 0.5f to get the area.
    */
    public static boolean getWindiness(List<Line> path, boolean b) {
      Line lineA = path.get(0);
      float sum = 0.0f;
      Vertex v1 = lineA.getVertex1();
      Vertex v2 = lineA.getVertex2();
      
      for(int i=1;i<path.size();i++) {
        
        sum += (v2.x-v1.x) * (v2.y+v1.y);
        v1=v2;
        Line lineB = path.get(i);
        v2 = lineB.getVertex2();

      }
      System.out.println("( windiness sum : " + sum + " )");

      if(sum<0) return true;
      else return false;
      
    }

    //make a shell Inset from path
    //this screws up the Line objects in path
    //so send a deep copy to avoid issues.

    public static List<Vertex> makeInsetPath(List<Line> path, boolean negateOffset, float offset) {

      List<Vertex> insetPath = new ArrayList<>();

      float shellOffset = offset;
      if(negateOffset) shellOffset = -shellOffset;

      //merge colinear segments
      path = merge(path);

      if(path.size()<2) return insetPath;

      boolean windiness = getWindiness(path,true); 
      if(!windiness) {
        Collections.reverse(path);

        for(Line l : path) {
          l.swapVerts();
        }
      }
      windiness = getWindiness(path,true);

     
      Line g = path.get(0);
      Vertex n = g.getNormal();
      //move line A along normal
      Line r = Line.add(g,n,shellOffset); 

      for(int i=1;i<path.size();i++) {
        
        Line g1 = path.get(i);
        Vertex n1 = g1.getNormal();
        //move line B along normal
        Line r1 = Line.add(g1,n1,shellOffset);

        //System.out.println("inset loop");
        //System.out.println(r + "  \n");
        //System.out.println(g1 + "::: normal => " + n1 + " = \n" + r1 +"\n");

        Vertex intersection = Line.getIntersection(r,r1,true);

        //System.out.println("intersection " + intersection);
        //System.out.println("");

        //Vertex intersection = g.getVertex1();

        insetPath.add(intersection);
        r = r1;
         
      }
      
      return insetPath;
    }
    public static void addGcodePrefix() {
      StringBuilder sb = new StringBuilder();
      
      String s = String.format("( minX: %f )\n", minX);
      sb.append(s);
      s = String.format("( maxX: %f )\n", maxX);
      sb.append(s);
      s = String.format("( minY: %f )\n", minY);
      sb.append(s);
      s = String.format("( maxY: %f )\n", maxY);
      sb.append(s);
      s = String.format("( minZ: %f )\n", minZ);
      sb.append(s);
      s = String.format("( maxZ: %f )\n", maxZ);
      sb.append(s);
      sb.append("G21\n");
      System.out.println(sb);
    }
    public static void convertToGcode(List<Line> pathIn, float z, float scale) {

      List<Line> path = new ArrayList<>(pathIn);
      StringBuilder sb = new StringBuilder();
      String fmt0 = "G0 X%.4f Y%.4f\n";
      String fmt = "G1 X%.4f Y%.4f\n";
      // move to first
      float x0 =  path.get(0).getVertex1().x * scale;
      float y0 =  path.get(0).getVertex1().y * scale;
      
      String s = String.format("G0 Z%.4f\n",z);
      sb.append(s);

      s = String.format(fmt0,x0,y0);
      sb.append(s);
      path.remove(0);

      for(Line l : path) {
         s = String.format("G1 X%.4f Y%.4f\n", l.getVertex1().x*scale ,l.getVertex1().y * scale);
         sb.append(s);
      }
      int end = path.size()-1;
      float x1 = path.get(end).getVertex2().x * scale;
      float y1 = path.get(end).getVertex2().y * scale;
      s = String.format(fmt, x1, y1);
      sb.append(s);
      s = String.format(fmt, x0, y0);
      sb.append(s);

      System.out.println(sb);

    }
    public static void convertVertexToGcode(List<Vertex> pathIn, float z) {

      List<Vertex> path = new ArrayList<>(pathIn);
      StringBuilder sb = new StringBuilder();
      String fmt0 = "G0 X%.4f Y%.4f\n";
      String fmt = "G1 X%.4f Y%.4f\n";

      // move to first
      if(pathIn.size()==0) return;

      float x0 =  path.get(0).x;
      float y0 =  path.get(0).y;
      
      String s = String.format("G0 Z%.4f\n",z);
      sb.append(s);

      s = String.format(fmt0,x0,y0);
      sb.append(s);
      path.remove(0);

      for(Vertex l : path) {
         s = String.format("G1 X%.4f Y%.4f\n",l.x ,l.y);
         sb.append(s);
      }
   
      s = String.format(fmt, x0, y0);
      sb.append(s);

      System.out.println(sb);

    }
    public static void fillGCode(List<Line> fill) {

      String fmt0 = "G0 X%.4f Y%.4f\n";
      String fmt1 = "G1 X%.4f Y%.4f\n";

      StringBuilder sb = new StringBuilder();

      if(fill.size() < 1) return;

      float x0 = fill.get(0).v1.x;
      float y0 = fill.get(0).v1.y;
      float x1 = fill.get(0).v2.x;
      float y1 = fill.get(0).v2.y;

      sb.append(String.format(fmt0, x0, y0));
      sb.append(String.format(fmt1, x1, y1));
      fill.remove(0);
      boolean reverse=false;
      while(fill.size()>0) {
          x0 = fill.get(0).v1.x;
          y0 = fill.get(0).v1.y;
          x1 = fill.get(0).v2.x;
          y1 = fill.get(0).v2.y;
         if(reverse) {
           sb.append(String.format(fmt1, x1, y1));
           sb.append(String.format(fmt1, x0, y0));
         }
         else {
           sb.append(String.format(fmt1, x0, y0));
           sb.append(String.format(fmt1, x1, y1));
         }
         reverse = !reverse;
         fill.remove(0);
      }
      System.out.println(sb);

    }
    public static void addFillX(List<List<Line>> paths, boolean writeOutput) {
      List<Line> fillList0 = new ArrayList<Line>();
      List<Line> fillList1 = new ArrayList<Line>();

      float dxl = 0.4f;
      for(float xl=minX; xl<=maxX; xl+=dxl) {
        
         Vertex a = new Vertex(xl,-200,0);
         Vertex b = new Vertex(xl,200,0);
         Line fillLine = new Line(a,b);
         List<Vertex> verts = new ArrayList<>();
         // find all intersects with the paths(if any)
         int pathNo = 0;
         for(List<Line> fPath : paths) {

            for (Line g : fPath) {

              Vertex v = g.findIntersection(fillLine);
              
              if(v!=null) {
                 v.pathNo = pathNo;
                 verts.add(v);
              }

            } // one path done

            // found all intersects now sort for fill order
            pathNo++;

         }//all paths considered

         //this only handles up to four intersections
         //TODO: fix
         if(verts.size()>1) {
              Collections.sort(verts,new SortByY());
              //System.out.println(verts);
              
              if(verts.size()==2 || verts.size()==4) {
                Line f0 = new Line(verts.get(0),verts.get(1));
                fillList0.add(f0);
                
                if(verts.size()==2) { 
                  ShellBuildDirection sbd = new ShellBuildDirection();
                  sbd.set(verts.get(0).pathNo, "INNER");
                  shell.add(sbd);
                }

              }
              if(verts.size()==4 || verts.size()==6) {
                Line f0 = new Line(verts.get(2),verts.get(3));
                fillList1.add(f0);

                ShellBuildDirection sbd0 = new ShellBuildDirection();
                sbd0.set(verts.get(0).pathNo, "INNER");
                shell.add(sbd0);

                if(verts.size()==4) { 
                  ShellBuildDirection sbd = new ShellBuildDirection();
                  sbd.set(verts.get(2).pathNo, "OUTER");
                  shell.add(sbd);
                }

              }
              
              if(verts.size()==6) {

                //Line f0 = new Line(verts.get(4),verts.get(5));
                //fillList2.add(f0);

                ShellBuildDirection sbd0 = new ShellBuildDirection();
                sbd0.set(verts.get(0).pathNo, "INNER");
                shell.add(sbd0);

                ShellBuildDirection sbd1 = new ShellBuildDirection();
                sbd1.set(verts.get(2).pathNo, "OUTER");
                shell.add(sbd1);

                if(verts.size()==6) {
                  ShellBuildDirection sbd = new ShellBuildDirection();
                  sbd.set(verts.get(4).pathNo, "OUTER");
                  shell.add(sbd);

                }
              }//six
         }
      }//all fill

      //play out fillList
      if(writeOutput) {
        fillGCode(fillList0);
        fillGCode(fillList1);
      }

    }
    public static List<Line> deepCopyPath(List<Line> path) {
      List<Line> r = new ArrayList<>();
      
      for(Line line : path) {
        r.add(new Line(line));
      }
      
      return r;
    }
    public static void main(String[] args) {

      Vertex v1 = new Vertex(-4,10,0);
      Vertex v2 = new Vertex(-4,0,0);
      Vertex v3 = new Vertex(-3,3.9,0);
      Vertex v4 = new Vertex(-4,4,0);

      v1 = new Vertex ( -11.190989,-1.5877931,0 );
      v2 = new Vertex ( -11.812888,-2.443747,0 );
      v3 = new Vertex ( -11.81289,-2.4437494,0 );
      v4 = new Vertex ( -11.88199,-2.5388556,0 );

      Line w = new Line(v1,v2);
      Line ww = new Line(v3,v4);
      Vertex vv = Line.getIntersection(ww,w,true);
      
      //Vertex n = w.getNormal();
      
      //System.out.println(vv);
      //w = Line.add(w,n);
      //System.out.println(w);
      // get vertex
      // then triangles
      
      List<Vertex> verts = new ArrayList<>();
      Set<Tri> tris = new HashSet<>();
      List<Line> lines = new ArrayList<>();
      List<Line> filteredLines = new ArrayList<>();

      Scanner scanner = new Scanner(System.in);
      
      int no_verts = scanner.nextInt();
      int no_tris = scanner.nextInt();
      scanner.nextLine();

      for(int i=0;i<no_verts;i++) {
        
        float x = scanner.nextFloat();
        float y = scanner.nextFloat();
        float z = scanner.nextFloat();
        scanner.nextLine();
        
        verts.add(new Vertex(x,y,z));

        if(x<minX) minX=x;
        if(x>maxX) maxX=x;

        if(y<minY) minY=y;
        if(y>maxY) maxY=y;

        if(z<minZ) minZ=z;
        if(z>maxZ) maxZ=z;
        
      }
      //System.out.println(verts);block2

      for(int i=0;i<no_tris;i++) {
        int s = scanner.nextInt();
        int vt1 = scanner.nextInt();
        int vt2 = scanner.nextInt();
        int vt3 = scanner.nextInt();
        scanner.nextLine();

        //System.out.println(verts.get(vt1));
        //System.out.println(verts.get(vt2));
        //System.out.println(verts.get(vt3));
        //System.out.println("");
        
        Tri tri = new Tri(verts.get(vt1),verts.get(vt2),verts.get(vt3));
        tris.add(tri);

        //STLwriter.writeFacet(tri);
      }

      float zoffset = minZ;
      float dz = 0.3f;
      float z = dz;
      

      addGcodePrefix();
      do {
        int i=0;
        for (Tri tr : tris) {
          //System.out.print("\n" + i + ": ");
          Line line = tr.getIntersects(z + zoffset);
          if(tr.intersectFound) {
              lines.add(line);

              //STLwriter.writeFacet(tr);

          }
          ++i; 
        }
      

      scanner.close();  
      for(Line l : lines) {
          if(l.toString().equals("Line: null to null")) {

          } else {
              filteredLines.add(l);
              //System.out.println(l);
          }
      }

      //System.out.println(filteredLines);
      //System.out.println(filteredLines.size());
      

      List<List<Line>> paths = new ArrayList<>();
      
      boolean newcode = true;

      if(newcode) {
      
        while(filteredLines.size()>2) {
          filteredLines = buildPath(filteredLines, z>3.5f && z<3.7f);
          paths.add(p);
          p = new ArrayList<>();
        }

        //if(paths.size() > 1)
        //    System.out.println(paths.get(1));
        //System.out.println(paths);

      }
      else {
        //seed the start
        p.add(filteredLines.get(0));
        filteredLines.remove(0);
        int s= filteredLines.size(); 
        for(int j=0;j<s;j++) { 
          filteredLines = findPath(filteredLines,p.size()-1,false);
        }
        paths.add(p);
        p = new ArrayList<>();
 
        if(filteredLines.size() > 1) {
          p.add(filteredLines.get(0));
          filteredLines.remove(0);
          //System.out.println(p);
          //System.out.println(filteredLines);
          int size = filteredLines.size();

          for(int j=0;j<size;j++) {
            filteredLines = findPath(filteredLines,p.size()-1, false);
          }
          paths.add(p);
        }
      
        //System.out.println(paths); 
      }//else 
      
      //fill
      //use fill algorithm to figure out shell inner or outer

      addFillX(paths, false);
      System.out.println("(" + shell + ")");

 
      System.out.println("( paths: " + paths.size() + ")");
      int pathNo = 0;
      boolean isOuter = false;
      for(List<Line> path : paths) {
        System.out.println("( PathNo: " + pathNo + " )");

        //remove isolated lines
        if(path.size() > 1) {

          convertToGcode(path,z,1.0f);
          System.out.println("(inset)");
          //get the paths inner or outer -ness  
          String Shellness = shell.findByPathNo(pathNo);
          if (Shellness.equals("NONE")) {
            isOuter=false;
            System.out.println("(shellness None for: " + pathNo + ")");
          }
          if (Shellness.equals("OUTER")) {
            isOuter = true;
          }
          else {
            isOuter = false;
          }
          List<Line> tempPath = deepCopyPath(path);
          convertVertexToGcode(makeInsetPath(tempPath,isOuter,0.4f),z);
          tempPath = deepCopyPath(path);
          convertVertexToGcode(makeInsetPath(tempPath,isOuter,0.8f),z);
          tempPath = deepCopyPath(path);
          convertVertexToGcode(makeInsetPath(tempPath,isOuter,1.2f),z);
          //convertToGcode(path,z,0.95f);
          //convertToGcode(path,z,0.90f);

        }
        pathNo++;
      }
      paths.clear();
      lines.clear();
      z+=dz;
      }while(z<(maxZ - minZ - dz));

      
    }//main
    


}
