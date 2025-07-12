#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        
            Bounds3 centroidBounds;
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds =
                    Union(centroidBounds, objects[i]->getBounds().Centroid());
            //get a big cube--centroidBounds
            int dim = centroidBounds.maxExtent();
            switch (dim) {
            case 0:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }
        if(splitMethod == SplitMethod::NAIVE){

            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);
        }
        else if(splitMethod == SplitMethod::SAH){
            int buckt = 12;
            std::vector<std::vector<Object*>> buckets(buckt);
            for (int i = 0; i < objects.size(); ++i) {
                //offset: get the position of now object in the big cube
                // and then, get its's bucket
                int b = buckt * centroidBounds.Offset(objects[i]->getBounds().Centroid())[dim];
                if (b >= buckt) b = buckt - 1;
                if(b < 0) b=0;
                buckets[b].push_back(objects[i]);
            }
            //now has put all objects in buckets
            int count = 0;
            int minCost = std::numeric_limits<float>::max();
            int minCostSplitBucket = 0;
            //[o,i] (i,buckt-1]
            auto UnionBuckt = [](std::vector<Object*>& VOB){
                Bounds3 bounds;
                for (int i = 0; i < VOB.size(); ++i) {
                    bounds = Union(bounds, VOB[i]->getBounds());
                }
                return bounds;
            };
            for (int i = 0; i < buckt - 1; ++i) {
                int count0 = 0, count1 = 0;
                Bounds3 bounds0, bounds1;
                for (int j = 0; j <= i; ++j) {
                    count0 += buckets[j].size();
                    bounds0 = Union(bounds0, UnionBuckt(buckets[j]));
                }
                for (int j = i + 1; j < buckt; ++j) {
                    count1 += buckets[j].size();
                    bounds1 = Union(bounds1, UnionBuckt(buckets[j]));
                }
                float SA=bounds0.SurfaceArea();
                float SB=bounds1.SurfaceArea();
                float cost = count0 * SA + count1 * SB+0.125;
                if (cost < minCost) {
                    minCost = cost;
                    minCostSplitBucket = i;
                }
            }
            auto beginning = objects.begin();
            auto middling = objects.begin() + objects.size()*minCostSplitBucket/buckt;
            auto ending = objects.end();

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);
        }
    }
    return node;

}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    Intersection isect;
    const std::array<int, 3> dirIsNeg{ int(ray.direction.x > 0), int(ray.direction.y > 0), int(ray.direction.z > 0) };
    if (node == nullptr || !node->bounds.IntersectP(ray, ray.direction_inv,dirIsNeg))
        return isect;
    if (node->left == nullptr && node->right == nullptr) {
        return node->object->getIntersection(ray);
    }
    Intersection hit1 = getIntersection(node->left, ray);
    Intersection hit2 = getIntersection(node->right, ray);
    return hit1.distance < hit2.distance ? hit1 : hit2;

}