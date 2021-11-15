#pragma once


enum SearchType {
    UnknownSearchMethod,
    DefaultSearchMethod,
    LinearSearchMethod,
    BinarySearchMethod,
    InterpolationSearchMethod,
    ExponentialSearchMethod,
    LinearRegressionSearchMethod
};

inline std::string toString(SearchType s) {
 std::string newMethod;
 switch(s) {
     case SearchType::LinearSearchMethod:
         newMethod = "LinearSearch";
             break;
     case SearchType::BinarySearchMethod:
             newMethod = "BinarySearch";
             break;
     case SearchType::InterpolationSearchMethod:
             newMethod = "InterpolationSearch";
             break;
     case SearchType::ExponentialSearchMethod:
         newMethod = "ExponentialSearch";
             break;
     case SearchType::DefaultSearchMethod:
         newMethod = "DefaultSearch";
         break;
     case SearchType::LinearRegressionSearchMethod:
         newMethod = "LinearRegressionSearch";
         break;
     default:
         break;
 }
 return newMethod;
}

enum PartitionType {UnknownPartition, SoAPartition, TreePartition, HashPartition};

inline std::string toString(PartitionType ptype) {
    switch (ptype) {
        case PartitionType::SoAPartition:
            return "SoA";
        case PartitionType::TreePartition:
            return "Tree";
        case PartitionType::HashPartition:
            return "Hash";
        default:
            return "Unknown";
    }
}

//enum Distribution { uniform, normal, zipf, wiki, exponential, books };
enum Distribution { uniform_dense, uniform_sparse, wiki, books, osm };

inline std::string toString(const Distribution d) {
    switch(d) {
    case Distribution::uniform_dense:
        return "uniform_dense";
    case Distribution::uniform_sparse:
        return "uniform_sparse";
    case Distribution::wiki:
        return "wiki";
    case Distribution::books:
        return "books";
    case Distribution::osm:
        return "osm";
    default:
        throw std::invalid_argument("Unsupported distribution.");
    }
}

inline Distribution distributionTypeFromString(const std::string& s) {
    if (s == "uniform_dense")
        return Distribution::uniform_dense;
    if (s == "uniform_sparse")
        return Distribution::uniform_sparse;
    if (s == "wiki")
        return Distribution::wiki;
    if (s == "books")
        return Distribution::books;
    if (s == "osm")
        return Distribution::osm;
    else
        throw std::invalid_argument("Unsupported distribution.");
}


enum Workload_kind {
    Point,
    Range,
    Insert
};


enum MutationType {
    PartitionTypeChangeMutationType,
    PartitionSearchChangeMutationType,
    EntrySearchChangeMutationType,
    HorizontalMergeMutationType,
    VerticalMergeMutationType,
    HorizontalSplitMutationType,
    VerticalSplitMutationType,
    KeySortedArraySoAPartitionMutationType,
    KeyHashPartitionMutationType,
    KeyTreePartitionMutationType,
    PartitionDeletionKeepRandomChildMutationType,
    PartitionDeletionRotateLeftMutationType,
    RotateLeftMutationType,
    ChainMutationType,
    LocalChainMutationType
};

inline std::string toString(MutationType t) {
    switch(t) {
        case MutationType::PartitionTypeChangeMutationType:
            return "PartitionTypeChangeMutation";
        case MutationType::PartitionSearchChangeMutationType:
            return "PartitionSearchChangeMutation";
        case MutationType::EntrySearchChangeMutationType:
            return "EntrySearchChangeMutation";
        case MutationType::HorizontalMergeMutationType:
            return "HorizontalMergeMutation";
        case MutationType::VerticalMergeMutationType:
            return "VerticalMergeMutation";
        case MutationType::HorizontalSplitMutationType:
            return "HorizontalSplitMutation";
        case MutationType::VerticalSplitMutationType:
            return "VerticalSplitMutation";
        case MutationType::KeySortedArraySoAPartitionMutationType:
            return "InsertSoAPartitionMutation";
        case MutationType::KeyHashPartitionMutationType:
            return "InsertHashPartitionMutation";
        case MutationType::KeyTreePartitionMutationType:
            return "InsertTreePartitionMutation";
        case MutationType::PartitionDeletionKeepRandomChildMutationType:
            return "PartitionDeletionKeepRandomChildMutation";
        case MutationType::PartitionDeletionRotateLeftMutationType:
            return "PartitionDeletionRotateLeftMutation";
        case MutationType::RotateLeftMutationType:
            return "RotateLeftMutation";
        case MutationType::ChainMutationType:
            return "ChainMutation";
        case MutationType::LocalChainMutationType:
            return "LocalChainMutation";
        default:
            return "UnknownMutationType";
    }
}
