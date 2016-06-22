void createSPARSOuterLoop();

void SparseCriteria::createSPARSOuterLoop()
{
  std::size_t indent = 2;

  // Reset parameters
  setup();
  visualizeOverlayNodes_ = false;  // DO NOT visualize all added nodes in a separate window
  denseCache_->resetCounters();

  // Get the ordering to insert vertices
  std::list<WeightedVertex> vertexInsertionOrder;
  getVertexInsertionOrdering(vertexInsertionOrder);

  // Error check order creation
  assert(vertexInsertionOrder.size() == getNumVertices() - queryVertices_.size());

  // Attempt to insert the vertices multiple times until no more succesful insertions occur
  secondSparseInsertionAttempt_ = false;
  std::size_t loopAttempt = 0;
  std::size_t sucessfulInsertions = 1;  // start with one so that while loop works
  while (sucessfulInsertions > 0)
  {
    std::cout << "Attempting to insert " << vertexInsertionOrder.size() << " vertices for the " << loopAttempt
              << " loop" << std::endl;

    // Sanity check
    if (loopAttempt > 3)
      OMPL_WARN("Suprising number of loop when attempting to insert nodes into SPARS graph: %u", loopAttempt);

    // Benchmark runtime
    time::point startTime = time::now();

    // ----------------------------------------------------------------------
    // Attempt to insert each vertex using the first 3 criteria
    if (!createSPARSInnerLoop(vertexInsertionOrder, sucessfulInsertions))
      break;

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);

    // Visualize
    if (visualizeSparseGraph_)
      visual_->viz1()->trigger();

    std::cout << "Succeeded in inserting " << sucessfulInsertions << " vertices on the " << loopAttempt
              << " loop, remaining uninserted vertices: " << vertexInsertionOrder.size()
              << " loop runtime: " << duration << " sec" << std::endl;
    loopAttempt++;

    // Increase the sparse delta a bit, but only after the first loop
    if (loopAttempt == 1)
    {
      // sparseDelta_ = getSecondarySparseDelta();
      std::cout << std::string(indent + 2, ' ') << "sparseDelta_ is now " << sparseDelta_ << std::endl;
      secondSparseInsertionAttempt_ = true;

      // Save collision cache, just in case there is a bug
      denseCache_->save();
    }

    bool debugOverRideJustTwice = true;
    if (debugOverRideJustTwice && loopAttempt == 1)
    {
      OMPL_WARN("Only attempting to add nodes twice for speed");
      break;
    }
  }

  // If we haven't animated the creation, just show it all at once
  if (!visualizeSparseGraph_)
  {
    sg_->displayDatabase(true, indent+4);
  }
  else if (sg_->visualizeSparseGraphSpeed_ < std::numeric_limits<double>::epsilon())
  {
    visual_->viz1()->trigger();
    usleep(0.001 * 1000000);
  }
}
