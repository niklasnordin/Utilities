//
//  ethanolCalculatorAppDelegate.h
//  ethanolCalculator
//
//  Created by Niklas Nordin on 2011-03-13.
//  Copyright 2011 nequam. All rights reserved.
//

#import <UIKit/UIKit.h>

@class RootViewController;

@class DetailViewController;

@interface ethanolCalculatorAppDelegate : NSObject <UIApplicationDelegate> {

}

@property (nonatomic, retain) IBOutlet UIWindow *window;

@property (nonatomic, retain) IBOutlet UISplitViewController *splitViewController;

@property (nonatomic, retain) IBOutlet RootViewController *rootViewController;

@property (nonatomic, retain) IBOutlet DetailViewController *detailViewController;

@end
