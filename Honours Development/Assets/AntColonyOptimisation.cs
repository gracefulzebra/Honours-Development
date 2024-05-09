using System;
using System.Collections.Generic;
using UnityEngine;

public class Node
{
    public Node() { }
    public Node(int _id) { id = _id; }

    public void AddNeighbourNode(int _neighbourIndex)
    {
        neighbourNodes.Add(_neighbourIndex);
    }

    public int id;
    public List<int> neighbourNodes = new List<int>();
}

public class Ant
{
    public Ant(int startingPos) { pos = startingPos; visitedNodes.Add(startingPos); }
    public int pos;
    public bool alive = true;
    public bool stuck = false;
    public bool ended = false;
    public List<int> visitedNodes = new List<int>();
    public float pathLength = 0f;
}

public class AntColonyOptimisation : MonoBehaviour
{
    [Header("Debug")]
    [SerializeField] private bool visibilityRadiusDebug = false;
    [SerializeField] private bool movementLineDebug = false;


    int numCities = 30;
    float[,] distanceMatrix; // = new float[numCities, numCities];
    float[,] pheromoneMatrix;
    Node[] nodes;
    [Range(0,100)]
    [SerializeField] private int numAnts = 1;
    Ant[] ants;

    float alpha = 0.2f;
    float beta = 0.8f;
    float rho = 0.2f;
    int maxIterations = 100; //Iterations relate to individual ants total steps
    int currentIteration = 0;
    int maxGeneration = 30; //Generations relate to each time the ants have completed the journey and new ants are to be made.
    int currentGeneration = 0;
    float bestDistance = float.MaxValue;

    bool solved = false;

    [Range(0, 100)]
    public float visibilityRadius;
    public int startingIndex = 1;
    public int finalTargetIndex = 11;
    Vector3[] debugLineData;

    [Range(2.0f, 10.0f)]
    public float tickRate = 2.0f;
    float time = 0.0f;
    bool started = false;
    
    public List<Transform> positionData = new List<Transform>();

    int status = 1;

    private void Awake()
    {
        InitialiseDatastructures();
        
    }

    private void Update()
    {
        if(Input.GetKeyDown(KeyCode.Space)) 
        { 
            started = !started; 
            if(started) { print("Started Algorithm"); }
            if(!started) { print("Stopped Algorithm"); }
        }

        if(started) { time += Time.deltaTime; }
        if(time > tickRate)
        {
            time = 0.0f;
            print("tick");
            InitialiseDistances(); 
            RunAlgorithm();
        }
    }

    void RunAlgorithm()
    {
        AlgorithmLoop();
        if(solved == false) { RunAlgorithm(); }
    }


    private void InitialiseDatastructures()
    {
        numCities = positionData.Count;

        distanceMatrix = new float[numCities, numCities];
        PropagateDistances();

        pheromoneMatrix = new float[numCities, numCities];
        PropagatePheromones();

        nodes = new Node[numCities];
        PropagateNodes();
    }

    private void InitialiseDistances()
    {
        distanceMatrix = new float[numCities, numCities];
        PropagateDistances();
        nodes = new Node[numCities];
        PropagateNodes();
    }

    void PropagateDistances()
    {
        for (int i = 0; i < numCities; i++)
        {
            for(int j = 0; j < numCities; j++)
            {
                //If node is comparing itself or if distance already calculated as ij == ji.
                if (i == j) { continue; }
                if (distanceMatrix[j, i] != 0)
                { 
                    distanceMatrix[i, j] = distanceMatrix[j, i];
                }
                else
                {
                    distanceMatrix[i, j] = CalculateDistance(positionData[i].position, positionData[j].position);
                }

            }
        }
    }

    void PropagatePheromones()
    {
        for(int i = 0; i < numCities; i++)
        {
            for(int j = 0; j < numCities; j++)
            {
                if (i == j) { continue; }
                float randomValue = UnityEngine.Random.Range(0.001f, 0.1f);
                pheromoneMatrix[i, j] = randomValue;
            }
        }
    }

    void PropagateNodes()
    {
        for(int i = 0; i < numCities; i++)
        {
            Node newNode = new Node(i);
            
            for(int j = 0; j < numCities; j++)
            {
                if (i == j) { continue; }
                if (distanceMatrix[i,j] < visibilityRadius) { newNode.AddNeighbourNode(j); }
            }
            nodes[i] = newNode;
        }
    }

    float CalculateDistance(Vector3 pointA, Vector3 pointB)
    {
        return (pointB - pointA).magnitude;
    }

    void AlgorithmLoop()
    {
        bool finish = false;
        status = 1;
        while (!finish)
        {
            switch(status)
            {
                case 1:

                    
                    InitialiseAnts();
                    currentIteration = 0;
                    status = 2;

                    break;

                case 2:

                    CalculateProbability();
                    
                    currentIteration++;

                    if (currentIteration <= maxIterations) { status = 2; }
                    else { status = 3; }

                    break;

                case 3:

                    EvaporatePheromone();
                    
                    currentGeneration++;

                    if(currentGeneration <= maxGeneration) { status = 1; }
                    else { status = 0; }

                    break;

                case 0:

                    ProcessFinalPath();
                    finish = true;
                    break;
            }
        }
    }    

    void InitialiseAnts()
    {
        ants = new Ant[numAnts];
        for(int i = 0; i <ants.Length; i++)
        {
            ants[i] = new Ant(startingIndex);
        }
    }

    void CalculateProbability()
    {
        for (int currentAntIndex = 0; currentAntIndex < ants.Length; currentAntIndex++) // For each ant
        {
            if (ants[currentAntIndex].stuck || ants[currentAntIndex].ended) { continue; }

            int currentNodePos = ants[currentAntIndex].pos;
            
            List<int> rawNeighbours = nodes[currentNodePos].neighbourNodes;
            for(int i = 0; i  < ants[currentAntIndex].visitedNodes.Count; i++)
            {
                if (rawNeighbours.Contains(ants[currentAntIndex].visitedNodes[i])) // Prune visited neighbours
                {
                    rawNeighbours.Remove(ants[currentAntIndex].visitedNodes[i]);
                }
            }

            int[] nodeNeighbourIndexes = rawNeighbours.ToArray();
            float[] nodeNeighbourProb = new float[nodeNeighbourIndexes.Length];

            if(nodeNeighbourIndexes.Length == 0) { ants[currentAntIndex].stuck = true; continue; }

            if(rawNeighbours.Contains(finalTargetIndex))  // Move to final node if available
            {
                ants[currentAntIndex].pos = finalTargetIndex;
                ants[currentAntIndex].visitedNodes.Add(finalTargetIndex);
                ants[currentAntIndex].ended = true;
                ants[currentAntIndex].pathLength += distanceMatrix[currentNodePos, finalTargetIndex];

                continue; 
            }

            for (int i = 0; i < nodeNeighbourIndexes.Length; i++) // For each connected node
            {
                float currentEdgePheromone = pheromoneMatrix[currentNodePos, nodeNeighbourIndexes[i]];
                float currentEdgeHeuristic = distanceMatrix[currentNodePos, nodeNeighbourIndexes[i]];

                float eqTop = MathF.Pow(currentEdgePheromone, alpha) * MathF.Pow(currentEdgeHeuristic, beta);
                float eqBot = 0;

                for (int j = 0; j < nodeNeighbourIndexes.Length; j++) // For all connected nodes
                {
                    int neighbourIndex = nodeNeighbourIndexes[j];
                    float edgePheromone = pheromoneMatrix[currentNodePos, neighbourIndex];
                    float edgeHeuristic = distanceMatrix[currentNodePos, neighbourIndex];

                    eqBot += MathF.Pow(edgePheromone, alpha) * MathF.Pow(edgeHeuristic, beta);
                }

                float totalProb = eqTop / eqBot;

                nodeNeighbourProb[i] = totalProb;
            }

            float[] cumulativeProb = new float[nodeNeighbourProb.Length];
            cumulativeProb[0] = nodeNeighbourProb[0];

            for (int k = 1; k < nodeNeighbourProb.Length; k++)
            {
                cumulativeProb[k] = cumulativeProb[k - 1] + nodeNeighbourProb[k];
            }

            float randomNumber = UnityEngine.Random.Range(0.0f, 1.0f);

            for(int m = 0; m < cumulativeProb.Length; m++) //Find next move for ant
            {
                if(randomNumber <= cumulativeProb[m])
                {
                    ants[currentAntIndex].pos = nodeNeighbourIndexes[m]; //Ant takes step to probabilistically chosen node
                    ants[currentAntIndex].visitedNodes.Add(nodeNeighbourIndexes[m]);
                    ants[currentAntIndex].pathLength += distanceMatrix[currentNodePos, nodeNeighbourIndexes[m]];
                    break;
                }
            }
        }
    }

    void EvaporatePheromone()
    {
        float[,] newPheromonesMatrix = new float[numCities, numCities];
        Array.Copy(pheromoneMatrix, newPheromonesMatrix, pheromoneMatrix.Length);

        for (int i = 0; i < numCities; i++)
        {
            for(int j = 0; j < numCities; j++)
            {
                if (i == j) { continue; }
                newPheromonesMatrix[i, j] = (1 - rho) * newPheromonesMatrix[i, j];
            }
        }

        for (int i = 0; i < numAnts; i++)
        {
            Ant currentAnt = ants[i];

            if (currentAnt.stuck) { currentAnt.pathLength = 100000; }
            else if(currentAnt.ended) { currentAnt.pathLength *= 0.5f; }

            int[] currentAntPath = currentAnt.visitedNodes.ToArray();

            for(int j = 0; j < currentAntPath.Length - 1; j++) 
            {
                newPheromonesMatrix[currentAntPath[j], currentAntPath[j + 1]] += 1 / currentAnt.pathLength;
            }
        }

        Array.Copy(newPheromonesMatrix, pheromoneMatrix, newPheromonesMatrix.Length);
    }

    void ProcessFinalPath()
    {
        List<Vector3> finalPath = new List<Vector3>();
        finalPath.Add(positionData[startingIndex].position);

        List<int> visitedNodes = new List<int>();

        int currentNodeIndex = 1;
        bool running = true;

        while (running)
        {
            int[] neighbourIndexes = nodes[currentNodeIndex].neighbourNodes.ToArray();

            if(neighbourIndexes.Length == 0) { solved = false; break; }

            float currentMax = 0f;
            int currentMaxNodeIndex = -1;
            for(int i = 0; i < neighbourIndexes.Length; i++)
            {
                if (visitedNodes.Contains(neighbourIndexes[i])) { continue; }
                if (pheromoneMatrix[currentNodeIndex, neighbourIndexes[i]] > currentMax)
                {
                    currentMax = pheromoneMatrix[currentNodeIndex, neighbourIndexes[i]];
                    currentMaxNodeIndex = neighbourIndexes[i];
                }
            }
            
            if(currentMaxNodeIndex == -1) { solved = false; break; }
            if(currentNodeIndex == finalTargetIndex) { finalPath.Add(positionData[finalTargetIndex].position); DebugAddPath(finalPath, visitedNodes); solved = true; break; }

            currentNodeIndex = currentMaxNodeIndex;

            finalPath.Add(positionData[currentNodeIndex].position);
            visitedNodes.Add(currentNodeIndex);
        }
    }

    void DebugAddPath(List<Vector3> debugPositionData, List<int> pathNodeIndex)
    {
        float pathLength = 0f;

        for (int i = 0; i < pathNodeIndex.Count - 1; i++)
        {
            pathLength += distanceMatrix[pathNodeIndex[i], pathNodeIndex[i + 1]];
        }

        if (pathLength < bestDistance)
        {
            bestDistance = pathLength;
            print("Current path length = " + bestDistance);
            debugLineData = debugPositionData.ToArray();
        }
    }

    private void OnDrawGizmos()
    {
        if(visibilityRadiusDebug)
        {
            Gizmos.color = Color.white;
            Gizmos.DrawWireSphere(positionData[startingIndex].position, visibilityRadius);
        }

        Gizmos.color = Color.red;
        Gizmos.DrawWireSphere(positionData[finalTargetIndex].position, 1.0f);

        if (movementLineDebug)
        {
            Gizmos.color = Color.blue;
            Gizmos.DrawLineStrip(debugLineData, false);
        }
    }

    
}
